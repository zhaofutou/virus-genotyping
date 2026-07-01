/**
 * Unified Virus Genotyping Server
 * Combines online (RIVM) + local (reference-based) subtyping.
 */

const fs = require('fs');
const path = require('path');
const express = require('express');
const multer = require('multer');
const { v4: uuidv4 } = require('uuid');
const { runPipeline } = require('./lib/pipeline');

const app = express();
const PORT = process.env.PORT || 3000;
const DATA_ROOT = path.join(__dirname, 'data');

// Middleware
app.use((req, res, next) => {
  console.log(`[${new Date().toISOString()}] ${req.method} ${req.url}`);
  next();
});
app.use(express.static(path.join(__dirname, 'public')));
app.use(express.json());

// Multer upload config
const upload = multer({
  storage: multer.diskStorage({
    destination: (req, file, cb) => {
      const dir = path.join(DATA_ROOT, req.jobId);
      fs.mkdirSync(dir, { recursive: true });
      cb(null, dir);
    },
    filename: (req, file, cb) => {
      if (file.originalname.endsWith('.zip')) {
        cb(null, 'upload.zip');
      } else {
        cb(null, file.originalname);
      }
    },
  }),
  limits: { fileSize: 50 * 1024 * 1024, files: 500 },
  fileFilter: (req, file, cb) => {
    const ext = path.extname(file.originalname).toLowerCase();
    if (['.seq', '.fasta', '.fa', '.zip'].includes(ext)) {
      cb(null, true);
    } else {
      cb(new Error(`Unsupported file type: ${ext}`));
    }
  },
});

// Job state
const jobs = new Map();

// Upload endpoint
app.post('/api/upload', (req, res, next) => {
  const jobId = uuidv4();
  req.jobId = jobId;
  next();
}, upload.any(), async (req, res) => {
  const jobId = req.jobId;
  const jobDir = path.join(DATA_ROOT, jobId);

  // Move non-zip files to sequences subdir
  if (req.files && req.files.length > 0) {
    const seqDir = path.join(jobDir, 'sequences');
    fs.mkdirSync(seqDir, { recursive: true });
    for (const f of req.files) {
      if (!f.originalname.endsWith('.zip') && path.dirname(f.path) !== seqDir) {
        fs.renameSync(f.path, path.join(seqDir, path.basename(f.path)));
      }
    }
  }

  const job = {
    status: 'queued',
    stage: 'starting',
    message: 'Job created',
    zipPath: null,
    clients: new Set(),
    result: null,
  };
  jobs.set(jobId, job);

  res.json({
    jobId,
    statusUrl: `/api/status/${jobId}`,
    downloadUrl: `/api/download/${jobId}`,
  });

  // Run pipeline asynchronously
  try {
    const emit = (stage, message) => {
      job.stage = stage;
      job.message = message;
      if (stage === 'done') job.status = 'done';
      for (const client of job.clients) {
        client.write(`data: ${JSON.stringify({ stage, message, status: job.status })}\n\n`);
      }
    };

    const result = await runPipeline(jobDir, emit);
    job.zipPath = result.zipPath;
    job.result = result;
    job.status = 'done';

    const finalMsg = {
      stage: 'done',
      message: 'Results ready',
      status: 'done',
      downloadUrl: `/api/download/${jobId}`,
      result: {
        hasConventionNames: result.hasConventionNames,
        comparison: result.comparison,
        online: result.online,
        local: result.local.map(r => ({
          name: r.name,
          length: r.length,
          noro: r.noro,
          rota: r.rota,
        })),
      },
    };
    for (const client of job.clients) {
      client.write(`data: ${JSON.stringify(finalMsg)}\n\n`);
    }
  } catch (err) {
    console.error(`[${jobId}] Pipeline error:`, err);
    job.status = 'error';
    job.message = err.message;
    for (const client of job.clients) {
      client.write(`data: ${JSON.stringify({ stage: 'error', message: err.message })}\n\n`);
    }
  }
});

// Status SSE endpoint
app.get('/api/status/:jobId', (req, res) => {
  const job = jobs.get(req.params.jobId);
  if (!job) return res.status(404).json({ error: 'Job not found' });

  res.writeHead(200, {
    'Content-Type': 'text/event-stream',
    'Cache-Control': 'no-cache',
    Connection: 'keep-alive',
  });

  // Send current state immediately
  res.write(`data: ${JSON.stringify({
    stage: job.stage,
    message: job.message,
    status: job.status,
    downloadUrl: job.status === 'done' ? `/api/download/${req.params.jobId}` : null,
    result: job.result ? {
      hasConventionNames: job.result.hasConventionNames,
      comparison: job.result.comparison,
      online: job.result.online,
      local: job.result.local?.map(r => ({
        name: r.name, length: r.length, noro: r.noro, rota: r.rota,
      })),
    } : null,
  })}\n\n`);

  if (job.status === 'done' || job.status === 'error') {
    res.end();
    return;
  }

  job.clients.add(res);
  req.on('close', () => job.clients.delete(res));
});

// Download endpoint
app.get('/api/download/:jobId', (req, res) => {
  const job = jobs.get(req.params.jobId);
  if (!job || !job.zipPath) return res.status(404).json({ error: 'Results not found' });
  res.download(job.zipPath, 'typing_results.zip');
});

// Cleanup old jobs (1 hour)
setInterval(() => {
  const cutoff = Date.now() - 3600000;
  for (const [jobId, job] of jobs) {
    const dir = path.join(DATA_ROOT, jobId);
    if (fs.existsSync(dir)) {
      const stat = fs.statSync(dir);
      if (stat.birthtimeMs < cutoff) {
        fs.rmSync(dir, { recursive: true, force: true });
        jobs.delete(jobId);
      }
    }
  }
}, 600000);

// Start
fs.mkdirSync(DATA_ROOT, { recursive: true });
app.listen(PORT, () => {
  console.log(`\n  Virus Genotyping Server`);
  console.log(`  ──────────────────────────────────`);
  console.log(`  Online:  RIVM Norovirus + Rotavirus A`);
  console.log(`  Local:   Reference-based alignment`);
  console.log(`  Port:    ${PORT}`);
  console.log(`  URL:     http://localhost:${PORT}\n`);
});
