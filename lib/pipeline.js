/**
 * Unified pipeline: runs online + local typing in parallel, then summarizes.
 */

const fs = require('fs');
const path = require('path');
const AdmZip = require('adm-zip');
const { Worker } = require('worker_threads');
const os = require('os');

const { loadReferences, parseFasta, buildKmerIndex, generateNorovirusCSV, generateRotavirusCSV } = require('./local-typing');
const { runOnlineTyping } = require('./online-typing');
const { concatSequences } = require('./concat');
const { getUploadedIndices, summarizeOnline, summarizeLocal } = require('./summarize');

class WorkerPool {
  constructor(size, workerScript, workerData) {
    this.workers = [];
    this.queue = [];
    for (let i = 0; i < size; i++) {
      const worker = new Worker(workerScript, { workerData });
      const entry = { worker, isBusy: false, currentResolve: null, currentReject: null };
      this.workers.push(entry);

      worker.on('message', (msg) => {
        entry.isBusy = false;
        if (msg.status === 'success') entry.currentResolve(msg.result);
        else entry.currentReject(new Error(msg.error));
        entry.currentResolve = null;
        entry.currentReject = null;
        this._next();
      });
      worker.on('error', (err) => {
        if (entry.currentReject) entry.currentReject(err);
        entry.isBusy = false;
        entry.currentResolve = null;
        entry.currentReject = null;
        this._next();
      });
    }
  }

  runTask(task) {
    return new Promise((resolve, reject) => {
      this.queue.push({ task, resolve, reject });
      this._next();
    });
  }

  _next() {
    if (this.queue.length === 0) return;
    const free = this.workers.find(w => !w.isBusy);
    if (!free) return;
    const item = this.queue.shift();
    free.isBusy = true;
    free.currentResolve = item.resolve;
    free.currentReject = item.reject;
    free.worker.postMessage(item.task);
  }

  close() {
    for (const w of this.workers) w.worker.terminate();
  }
}

async function runPipeline(jobDir, emit) {
  const startedAt = Date.now();
  const elapsed = () => Math.round((Date.now() - startedAt) / 1000);

  try {
    // Step 0: Extract zip if present
    const zipFile = path.join(jobDir, 'upload.zip');
    const seqDir = path.join(jobDir, 'sequences');
    if (fs.existsSync(zipFile)) {
      emit('extracting', 'Extracting zip file...');
      new AdmZip(zipFile).extractAllTo(seqDir, true);
    } else {
      fs.mkdirSync(seqDir, { recursive: true });
    }

    // Step 1: Concatenate
    emit('concatenating', 'Concatenating sequences...');
    const { path: fastaPath, count } = concatSequences(seqDir, jobDir);
    emit('concatenating', `${count} sequences concatenated`);

    const fastaText = fs.readFileSync(fastaPath, 'utf8');
    const queries = parseFasta(fastaText);

    // Determine which indices were uploaded (for summary)
    const uploadedIndices = getUploadedIndices(seqDir);
    const hasConventionNames = uploadedIndices.size > 0;

    // Step 2: Run online + local typing in parallel
    emit('typing', `Running online (RIVM) + local typing for ${count} sequences...`);

    const [onlineResult, localResults] = await Promise.all([
      // Online typing (RIVM)
      runOnlineTyping(fastaText).catch(err => {
        emit('typing', `Online typing error: ${err.message}`);
        return { noro: { error: err.message }, rota: { error: err.message } };
      }),
      // Local typing
      (async () => {
        const projectRoot = path.resolve(__dirname, '..');
        emit('typing', 'Loading local references...');
        const allRefs = loadReferences(projectRoot);
        const noroRefs = allRefs.filter(r => r.virus === 'noro');
        const rotaRefs = allRefs.filter(r => r.virus === 'rota');
        const noroIdx = buildKmerIndex(noroRefs, 13);
        const rotaIdx = buildKmerIndex(rotaRefs, 13);

        const numCores = os.cpus().length;
        const poolSize = Math.max(1, numCores);
        const pool = new WorkerPool(poolSize, path.join(__dirname, 'worker.js'), { baseDir: projectRoot });

        let typedCount = 0;
        const options = {
          k: 13, candidateLimit: 40, minOverlap: 100,
          noroThresholds: { VP1: 80, RdRp: 85 },
          rotaThresholds: { VP7: 80, VP4: 85 },
        };

        const promises = queries.map(q => {
          return pool.runTask({ name: q.name, sequence: q.sequence, options })
            .then(result => {
              typedCount++;
              emit('typing', `Local typing: ${typedCount}/${queries.length} completed`);
              return result;
            });
        });

        const results = await Promise.all(promises);
        pool.close();
        return results;
      })(),
    ]);

    // Step 3: Generate CSVs
    emit('generating', 'Generating result files...');

    // Local CSVs
    const localNoroCSV = generateNorovirusCSV(localResults);
    const localRotaCSV = generateRotavirusCSV(localResults);
    fs.writeFileSync(path.join(jobDir, 'local_noro_result.csv'), localNoroCSV);
    fs.writeFileSync(path.join(jobDir, 'local_rota_result.csv'), localRotaCSV);

    // Online CSVs
    if (onlineResult.noro?.csvContent) {
      fs.writeFileSync(path.join(jobDir, 'online_noro_result.csv'), onlineResult.noro.csvContent);
    }
    if (onlineResult.rota?.csvContent) {
      fs.writeFileSync(path.join(jobDir, 'online_rota_result.csv'), onlineResult.rota.csvContent);
    }

    // Step 4: Summarize (only if naming convention is followed)
    let onlineSummary = null;
    let localSummary = null;

    if (hasConventionNames) {
      emit('summarizing', 'Building 1-72 summary tables...');

      onlineSummary = summarizeOnline(
        onlineResult.noro?.results || null,
        onlineResult.rota?.results || null,
        uploadedIndices,
      );
      localSummary = summarizeLocal(localResults, uploadedIndices);

      // Write summary CSVs
      const writeCSV = (filepath, data) => {
        let csv = 'index,typing_result\n';
        for (const r of data) csv += `${r.index},"${r.typing_result}"\n`;
        fs.writeFileSync(filepath, csv);
      };

      writeCSV(path.join(jobDir, 'online_noro_summary.csv'), onlineSummary.noroSummary);
      writeCSV(path.join(jobDir, 'online_rota_summary.csv'), onlineSummary.rotaSummary);
      writeCSV(path.join(jobDir, 'online_final_1_72.csv'), onlineSummary.master);

      writeCSV(path.join(jobDir, 'local_noro_summary.csv'), localSummary.noroSummary);
      writeCSV(path.join(jobDir, 'local_rota_summary.csv'), localSummary.rotaSummary);
      writeCSV(path.join(jobDir, 'local_final_1_72.csv'), localSummary.master);
    }

    // Step 5: Package
    emit('packaging', 'Packaging results...');
    const zip = new AdmZip();
    const resultFiles = [
      'concated2.fasta',
      'online_noro_result.csv', 'online_rota_result.csv',
      'local_noro_result.csv', 'local_rota_result.csv',
      'online_noro_summary.csv', 'online_rota_summary.csv', 'online_final_1_72.csv',
      'local_noro_summary.csv', 'local_rota_summary.csv', 'local_final_1_72.csv',
    ];
    for (const f of resultFiles) {
      const fp = path.join(jobDir, f);
      if (fs.existsSync(fp)) zip.addLocalFile(fp);
    }
    const zipPath = path.join(jobDir, 'results.zip');
    zip.writeZip(zipPath);

    emit('done', `Complete in ${elapsed()}s`);

    return {
      zipPath,
      onlineSummary: onlineSummary?.master || null,
      localSummary: localSummary?.master || null,
      hasConventionNames,
      online: {
        noro: onlineResult.noro?.results || [],
        rota: onlineResult.rota?.results || [],
        noroError: onlineResult.noro?.error || null,
        rotaError: onlineResult.rota?.error || null,
      },
      local: localResults,
    };
  } catch (err) {
    emit('error', err.message);
    throw err;
  }
}

module.exports = { runPipeline };
