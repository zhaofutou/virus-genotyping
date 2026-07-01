/**
 * NCBI BLAST-based typing.
 * Submits sequences to NCBI BLAST, parses top hits for genotyping.
 * When genogroup is known (GI/GII), uses targeted database for better results.
 */

const https = require('https');
const http = require('http');

const BLAST_API = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi';

/**
 * Submit a sequence to NCBI BLAST.
 * @param {string} sequence - nucleotide sequence (no header)
 * @param {object} opts
 * @param {string} [opts.program] - blastn (default)
 * @param {string} [opts.database] - nr/nt/refseq_rna (default: nr)
 * @param {string} [opts.genogroup] - GI, GII, etc. to narrow the search
 * @returns {Promise<string>} RID (Request ID)
 */
async function submitBlast(sequence, opts = {}) {
  const { program = 'blastn', database = 'nr', genogroup } = opts;

  // Build Entrez query to narrow results if genogroup is known
  let entrezQuery = '';
  if (genogroup) {
    // Search within Norovirus genogroup
    entrezQuery = `"Norovirus ${genogroup}"[Organism] OR "norovirus ${genogroup}"[Title]`;
  }

  const params = new URLSearchParams({
    CMD: 'Put',
    PROGRAM: program,
    DATABASE: database,
    QUERY: sequence,
    FORMAT_TYPE: 'JSON2',
    HITLIST_SIZE: '20',
    EXPECT: '10',
    WORD_SIZE: '11',
    MATRIX_NAME: 'BLOSUM62',
    GAPCOSTS: '5 2',
    COMPOSITION_BASED_STATISTICS: '2',
  });

  if (entrezQuery) {
    params.set('ENTREZ_QUERY', entrezQuery);
  }

  const body = params.toString();

  return new Promise((resolve, reject) => {
    const url = new URL(BLAST_API);
    const options = {
      hostname: url.hostname,
      path: url.pathname,
      method: 'POST',
      headers: {
        'Content-Type': 'application/x-www-form-urlencoded',
        'Content-Length': Buffer.byteLength(body),
        'User-Agent': 'virus-genotyping/1.0',
      },
    };

    const req = https.request(options, (res) => {
      let data = '';
      res.on('data', chunk => data += chunk);
      res.on('end', () => {
        // Parse RID from response
        const ridMatch = data.search(/RID\s*=\s*([A-Z0-9]+)/i) !== -1
          ? data.match(/RID\s*=\s*([A-Z0-9]+)/i)
          : null;
        // Try JSON first
        try {
          const json = JSON.parse(data);
          if (json.BlastOutput2?.[0]?.report?.results?.search?.RID) {
            return resolve(json.BlastOutput2[0].report.results.search.RID);
          }
        } catch {}

        if (ridMatch) {
          resolve(ridMatch[1]);
        } else {
          reject(new Error('BLAST submit failed: no RID returned. Response: ' + data.slice(0, 500)));
        }
      });
    });

    req.on('error', reject);
    req.setTimeout(30000, () => { req.destroy(); reject(new Error('BLAST submit timeout')); });
    req.write(body);
    req.end();
  });
}

/**
 * Poll BLAST status until results are ready.
 * @param {string} rid
 * @param {number} maxWait - max wait in ms (default 5 min)
 * @returns {Promise<boolean>}
 */
async function pollBlast(rid, maxWait = 300000) {
  const interval = 10000;
  const maxAttempts = Math.ceil(maxWait / interval);

  for (let i = 0; i < maxAttempts; i++) {
    await new Promise(r => setTimeout(r, interval));

    const status = await new Promise((resolve, reject) => {
      const url = `${BLAST_API}?CMD=Get&FORMAT_TYPE=JSON2&RID=${rid}`;
      https.get(url, { headers: { 'User-Agent': 'virus-genotyping/1.0' } }, (res) => {
        let data = '';
        res.on('data', chunk => data += chunk);
        res.on('end', () => {
          // Check if it's still running
          if (data.includes('Status=WAITING') || data.includes('Status=QUEUED') || data.includes('Status=RUNNING')) {
            resolve('running');
          } else if (data.includes('Status=FAILED')) {
            resolve('failed');
          } else if (data.includes('Status=UNKNOWN')) {
            resolve('unknown');
          } else if (data.includes('BlastOutput2') || data.includes('BlastOutput')) {
            resolve('ready');
          } else {
            resolve('unknown');
          }
        });
      }).on('error', reject);
    });

    if (status === 'ready') return true;
    if (status === 'failed') throw new Error(`BLAST job ${rid} failed`);
    if (status === 'unknown') throw new Error(`BLAST job ${rid} status unknown`);
    // Still running, continue polling
  }

  throw new Error(`BLAST job ${rid} timed out after ${maxWait / 1000}s`);
}

/**
 * Fetch and parse BLAST results.
 * @param {string} rid
 * @returns {Promise<object>} Parsed results
 */
async function fetchBlastResults(rid) {
  return new Promise((resolve, reject) => {
    const url = `${BLAST_API}?CMD=Get&FORMAT_TYPE=JSON2&RID=${rid}&ALIGNMENTS=20`;
    https.get(url, { headers: { 'User-Agent': 'virus-genotyping/1.0' } }, (res) => {
      let data = '';
      res.on('data', chunk => data += chunk);
      res.on('end', () => {
        try {
          const json = JSON.parse(data);
          const search = json.BlastOutput2?.[0]?.report?.results?.search;
          if (!search) return reject(new Error('No search results in BLAST output'));

          const hits = (search.hits || []).map(hit => {
            const desc = hit.description?.[0] || {};
            const hsps = hit.hsps?.[0] || {};
            return {
              title: desc.title || '',
              accession: desc.accession || '',
              taxid: desc.taxid || 0,
              sciname: desc.sciname || '',
              score: hsps.score || 0,
              evalue: hsps.evalue || 0,
              identity: hsps.identity || 0,
              alignLen: hsps.align_len || 0,
              percentIdentity: hsps.align_len ? ((hsps.identity / hsps.align_len) * 100).toFixed(2) : '0',
              queryFrom: hsps.query_from || 0,
              queryTo: hsps.query_to || 0,
              hitFrom: hsps.hit_from || 0,
              hitTo: hsps.hit_to || 0,
            };
          });

          resolve({
            queryLen: search.query_len || 0,
            rid,
            hits,
          });
        } catch (err) {
          reject(new Error('Failed to parse BLAST results: ' + err.message));
        }
      });
    }).on('error', reject);
  });
}

/**
 * Extract norovirus genotype from BLAST hit title.
 * Looks for patterns like "GI.3", "GII.4", "GI.P13", "GII.P16" in the title.
 */
function extractGenotype(title, sciname) {
  const text = `${title} ${sciname}`;

  // Match capsid genotype: GI.3, GII.4, GII.17, etc.
  const capsidMatch = text.match(/\bG(I{1,3}|IV|V{1,3}|IX|X)\.(\d+[a-z]?)\b/i);

  // Match polymerase genotype: GI.P13, GII.P16, etc.
  const polMatch = text.match(/\bG(I{1,3}|IV|V{1,3}|IX|X)\.P(\d+[a-z]?(?:\.\w+)?)\b/i);

  // Match combined: GI.3d or GII.4/Sydney
  const combinedMatch = text.match(/\bG(I{1,3}|IV|V{1,3}|IX|X)\.(\d+[a-z]?)(?:\/[\w]+)?\b/i);

  return {
    capsid: capsidMatch ? capsidMatch[0] : '',
    polymerase: polMatch ? polMatch[0] : '',
    combined: combinedMatch ? combinedMatch[0] : '',
    raw: text,
  };
}

/**
 * Extract rotavirus genotype from BLAST hit title.
 * Looks for G1-G12 and P[4]-P[19] patterns.
 */
function extractRotaGenotype(title, sciname) {
  const text = `${title} ${sciname}`;

  // VP7 G-type: G1, G2, ... G12
  const gMatch = text.match(/\bG(\d{1,2})\b/);

  // VP4 P-type: P[4], P[8], etc.
  const pMatch = text.match(/\bP\[(\d{1,2})\]/);

  return {
    gType: gMatch ? `G${gMatch[1]}` : '',
    pType: pMatch ? `P[${pMatch[1]}]` : '',
  };
}

/**
 * Run BLAST for a single sequence and extract genotyping info.
 * @param {string} name - sequence name
 * @param {string} sequence - nucleotide sequence
 * @param {object} opts
 * @param {string} [opts.genogroup] - known genogroup (GI/GII) for targeted search
 * @param {string} [opts.virus] - 'noro' or 'rota' (auto-detect if not set)
 * @returns {Promise<object>}
 */
async function blastType(name, sequence, opts = {}) {
  const { genogroup, virus } = opts;

  const rid = await submitBlast(sequence, { genogroup });
  await pollBlast(rid);
  const results = await fetchBlastResults(rid);

  // Extract genotypes from top hits
  const topHits = results.hits.slice(0, 10);
  let bestCapsid = '';
  let bestPol = '';
  let bestGType = '';
  let bestPType = '';

  for (const hit of topHits) {
    if (!bestCapsid || !bestPol) {
      const noro = extractGenotype(hit.title, hit.sciname);
      if (!bestCapsid && noro.capsid) bestCapsid = noro.capsid;
      if (!bestPol && noro.polymerase) bestPol = noro.polymerase;
    }
    if (!bestGType || !bestPType) {
      const rota = extractRotaGenotype(hit.title, hit.sciname);
      if (!bestGType && rota.gType) bestGType = rota.gType;
      if (!bestPType && rota.pType) bestPType = rota.pType;
    }
  }

  const topHit = results.hits[0];

  return {
    name,
    length: sequence.length,
    blast: {
      rid,
      topHitTitle: topHit?.title || '',
      topHitIdentity: topHit?.percentIdentity || '0',
      topHitEvalue: topHit?.evalue || 0,
      queryLen: results.queryLen,
      hitCount: results.hits.length,
    },
    noro: {
      capsid: bestCapsid,
      polymerase: bestPol,
    },
    rota: {
      gType: bestGType,
      pType: bestPType,
    },
  };
}

module.exports = { submitBlast, pollBlast, fetchBlastResults, blastType, extractGenotype, extractRotaGenotype };
