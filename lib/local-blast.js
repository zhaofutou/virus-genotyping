const fs = require('fs');
const path = require('path');
const { exec } = require('child_process');

const PROJECT_ROOT = '/home/ts_admin/virus-genotyping';
const DB_DIR = path.join(PROJECT_ROOT, 'data/blast_db');
const DB_PATH = path.join(DB_DIR, 'virus_db');

/**
 * Scan directory recursively for FASTA files
 */
function scanFastaFiles(dir) {
  const files = [];
  if (!fs.existsSync(dir)) return files;
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      files.push(...scanFastaFiles(full));
    } else {
      const ext = path.extname(entry.name).toLowerCase();
      if (['.fa', '.fasta'].includes(ext)) {
        files.push(full);
      }
    }
  }
  return files;
}

/**
 * Initialize and compile BLAST database if needed
 */
function initBlastDatabase() {
  if (!fs.existsSync(DB_DIR)) {
    fs.mkdirSync(DB_DIR, { recursive: true });
  }

  const dbFiles = ['.nhr', '.nin', '.nsq'];
  const dbExists = dbFiles.every(ext => fs.existsSync(DB_PATH + ext));

  // If DB exists, skip build unless references have changed
  if (dbExists) {
    console.log('[BLAST] Local database already exists.');
    return Promise.resolve();
  }

  console.log('[BLAST] Compiling local database from references...');
  const refDir = path.join(PROJECT_ROOT, 'reference');
  const refFiles = scanFastaFiles(refDir);
  if (refFiles.length === 0) {
    return Promise.reject(new Error('No reference sequences found under reference/'));
  }

  // Concatenate all reference files
  let concatenated = '';
  for (const file of refFiles) {
    concatenated += fs.readFileSync(file, 'utf8') + '\n';
  }

  const allRefsFasta = path.join(DB_DIR, 'all_references.fasta');
  fs.writeFileSync(allRefsFasta, concatenated);

  return new Promise((resolve, reject) => {
    const cmd = `makeblastdb -in "${allRefsFasta}" -dbtype nucl -out "${DB_PATH}"`;
    exec(cmd, (err, stdout, stderr) => {
      if (err) {
        console.error('[BLAST] makeblastdb error:', stderr);
        return reject(err);
      }
      console.log('[BLAST] makeblastdb success:', stdout.trim());
      resolve();
    });
  });
}

/**
 * Extract genotype from local BLAST hit title, using query name to filter by segment
 */
function parseGenotypeFromHit(title, queryName) {
  const text = title.toLowerCase();
  const queryLower = (queryName || '').toLowerCase();
  
  let virus = 'unknown';
  let genotype = '';

  if (text.includes('noro') || text.includes('nov-')) {
    virus = 'norovirus';
    const hasVp1 = queryLower.includes('vp1') || queryLower.includes('capsid');
    const hasRdrp = queryLower.includes('rdrp') || queryLower.includes('pol');
    
    // Match capsid (e.g., GII.4) or polymerase (e.g., GII.P16)
    const capsidMatch = title.match(/\bG(I{1,3}|IV|V{1,3}|IX|X)\.\d+[a-z]?\b/i);
    const polMatch = title.match(/\bG(I{1,3}|IV|V{1,3}|IX|X)\.P\d+[a-z]?\b/i);

    if (hasVp1 && capsidMatch) {
      genotype = capsidMatch[0];
    } else if (hasRdrp && polMatch) {
      genotype = polMatch[0];
    } else {
      const matches = title.match(/\bG(I{1,3}|IV|V{1,3}|IX|X)\.(P?\d+[a-z]?)\b/i);
      if (matches) genotype = matches[0];
    }
  } else if (text.includes('rota') || text.includes('rva-')) {
    virus = 'rotavirus';
    const hasV4 = queryLower.includes('v4') || queryLower.includes('vp4');
    const hasV7 = queryLower.includes('v7') || queryLower.includes('vp7');
    
    const gMatch = title.match(/G(\d{1,2})/i);
    const pMatch = title.match(/P\[?(\d{1,2})\]?/i);
    
    if (hasV4 && pMatch) {
      genotype = `P[${pMatch[1]}]`;
    } else if (hasV7 && gMatch) {
      genotype = gMatch[0];
    } else {
      if (gMatch && pMatch) genotype = `${gMatch[0]}P[${pMatch[1]}]`;
      else if (gMatch) genotype = gMatch[0];
      else if (pMatch) genotype = `P[${pMatch[1]}]`;
    }
  } else if (text.includes('sapo')) {
    virus = 'sapovirus';
    const matches = title.match(/GI{1,2}\.\d+/i);
    if (matches) genotype = matches[0];
  } else if (text.includes('astro')) {
    virus = 'astrovirus';
    const matches = title.match(/hastv-\d+/i);
    if (matches) genotype = matches[0].toUpperCase();
  } else if (text.includes('adeno')) {
    virus = 'adenovirus';
    const matches = title.match(/hadv-\d+/i);
    if (matches) genotype = matches[0].toUpperCase();
  }

  return { virus, genotype, rawTitle: title };
}

/**
 * Query sequence(s) against local BLAST database
 * @param {string} queryFastaPath
 * @returns {Promise<Array>} Results per query
 */
async function runLocalBlast(queryFastaPath) {
  await initBlastDatabase();

  return new Promise((resolve, reject) => {
    // -outfmt 15 outputs JSON format
    const cmd = `blastn -query "${queryFastaPath}" -db "${DB_PATH}" -outfmt 15 -max_target_seqs 5`;
    exec(cmd, { maxBuffer: 10 * 1024 * 1024 }, (err, stdout, stderr) => {
      if (err) {
        console.error('[BLAST] blastn error:', stderr);
        return reject(err);
      }

      try {
        const json = JSON.parse(stdout);
        const searchResults = json.BlastOutput2 || [];
        const parsedResults = [];

        for (const item of searchResults) {
          const search = item.report?.results?.search;
          if (!search) continue;

          const queryTitle = search.query_title || '';
          const queryLen = search.query_len || 0;
          const hits = search.hits || [];

          // Classify based on top hits
          let bestVirus = 'unknown';
          let bestGenotype = '';
          let topHitIdentity = '0';
          let topHitTitle = '';

          if (hits.length > 0) {
            const topHit = hits[0];
            topHitTitle = topHit.description?.[0]?.title || '';
            const hsps = topHit.hsps?.[0] || {};
            const alignLen = hsps.align_len || 0;
            topHitIdentity = alignLen ? ((hsps.identity / alignLen) * 100).toFixed(2) : '0';

            // Gather consensus genotype from first few hits
            for (const hit of hits.slice(0, 3)) {
              const hitTitle = hit.description?.[0]?.title || '';
              const parsed = parseGenotypeFromHit(hitTitle, queryTitle);
              if (parsed.genotype) {
                bestVirus = parsed.virus;
                bestGenotype = parsed.genotype;
                break;
              }
            }
          }

          parsedResults.push({
            name: queryTitle,
            length: queryLen,
            virus: bestVirus,
            genotype: bestGenotype,
            blast: {
              topHitTitle,
              topHitIdentity,
              hitCount: hits.length,
            }
          });
        }

        resolve(parsedResults);
      } catch (parseErr) {
        reject(new Error('Failed to parse local BLAST output: ' + parseErr.message));
      }
    });
  });
}

module.exports = { initBlastDatabase, runLocalBlast, parseGenotypeFromHit };
