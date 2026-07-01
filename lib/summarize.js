/**
 * Summarize typing results into a 1-72 index table.
 * Supports both online (RIVM CSV) and local (engine output) results.
 */

const path = require('path');

// Parse sample name: May05V4 → { month: 'May', index: 5, type: 'V4' }
// Types: V4=rotavirus VP4, V7=rotavirus VP7, GI=norovirus GI, GII=norovirus GII
function parseSampleName(filename) {
  const base = path.basename(filename).replace(/\.(seq|fasta|fa)$/i, '');
  const prefix = base.split('_(')[0];

  let type = null;
  if (prefix.includes('GII')) type = 'GII';
  else if (prefix.includes('GI') && !prefix.includes('GII')) type = 'GI';
  else if (prefix.includes('V4')) type = 'V4';
  else if (prefix.includes('V7')) type = 'V7';

  if (!type) return null;

  const cleanedPrefix = prefix.replace(new RegExp(type, 'g'), '');
  const matches = cleanedPrefix.match(/\d+/g);
  if (!matches) return null;

  for (const m of matches) {
    const idx = parseInt(m, 10);
    if (idx >= 1 && idx <= 72) {
      const monthMatch = prefix.match(/^[a-zA-Z]+/);
      const month = monthMatch ? monthMatch[0] : '';
      return { month, index: idx, type };
    }
  }
  return null;
}

// Determine which indices were uploaded and what types they have
function getUploadedIndices(seqDir) {
  const fs = require('fs');
  const indices = new Map(); // index -> Set<type>

  function scan(dir) {
    if (!fs.existsSync(dir)) return;
    for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
      const full = path.join(dir, entry.name);
      if (entry.isDirectory()) { scan(full); continue; }
      if (!/\.(seq|fasta|fa)$/i.test(entry.name)) continue;
      if (entry.name.startsWith('concated')) continue;

      const parsed = parseSampleName(entry.name);
      if (!parsed) continue;
      if (!indices.has(parsed.index)) indices.set(parsed.index, new Set());
      indices.get(parsed.index).add(parsed.type);
    }
  }
  scan(seqDir);
  return indices;
}

/**
 * Build 1-72 summary from online RIVM results.
 * @param {Array} noroResults - Rows from RIVM norovirus CSV/JSON
 * @param {Array} rotaResults - Rows from RIVM rotavirus CSV/JSON
 * @param {Map} uploadedIndices - index -> Set<type> from uploaded files
 */
function summarizeOnline(noroResults, rotaResults, uploadedIndices) {
  const noroBestCapsid = {};
  const noroBestPol = {};
  const rotaBestG = {};
  const rotaBestP = {};

  // Process norovirus online results
  if (noroResults) {
    for (const row of noroResults) {
      const name = row['Name'] || row['name'] || '';
      const parsed = parseSampleName(name);
      if (!parsed) continue;
      const idx = parsed.index;
      const score = parseFloat(row['BLAST score'] || row['blast_score'] || 0);
      const capsid = row['capsid type'] || row['capsid_1'] || '';
      const pol = row['polymerase type'] || row['polymerase_1'] || '';

      if (capsid && (!noroBestCapsid[idx] || score > noroBestCapsid[idx].score)) {
        noroBestCapsid[idx] = { genotype: capsid, score };
      }
      if (pol && (!noroBestPol[idx] || score > noroBestPol[idx].score)) {
        noroBestPol[idx] = { genotype: pol, score };
      }
    }
  }

  // Process rotavirus online results
  if (rotaResults) {
    for (const row of rotaResults) {
      const name = row['Name'] || row['name'] || '';
      const parsed = parseSampleName(name);
      if (!parsed) continue;
      const idx = parsed.index;
      const score = parseFloat(row['BLAST score'] || row['blast_score'] || 0);

      const val1 = row['cluster result'] || row['cluster/support_1'] || '';
      const val2 = row['cluster/support_2'] || '';

      let gGen = '', pGen = '';
      for (const val of [val1, val2]) {
        if (!val) continue;
        if (/^G\d+$/i.test(val)) gGen = val;
        if (/^P\[\d+\]$/i.test(val)) pGen = val;
      }

      if (parsed.type === 'V7' && gGen) {
        if (!rotaBestG[idx] || score > rotaBestG[idx].score) {
          rotaBestG[idx] = { genotype: gGen, score };
        }
      }
      if (parsed.type === 'V4' && pGen) {
        if (!rotaBestP[idx] || score > rotaBestP[idx].score) {
          rotaBestP[idx] = { genotype: pGen, score };
        }
      }
      if (parsed.type !== 'V4' && parsed.type !== 'V7') {
        if (gGen && (!rotaBestG[idx] || score > rotaBestG[idx].score)) {
          rotaBestG[idx] = { genotype: gGen, score };
        }
        if (pGen && (!rotaBestP[idx] || score > rotaBestP[idx].score)) {
          rotaBestP[idx] = { genotype: pGen, score };
        }
      }
    }
  }

  return buildMasterTable(uploadedIndices, {
    noroBestCapsid, noroBestPol, rotaBestG, rotaBestP,
  });
}

/**
 * Build 1-72 summary from local typing results.
 * @param {Array} localResults - Array of typeSequence() output objects
 * @param {Map} uploadedIndices - index -> Set<type>
 */
function summarizeLocal(localResults, uploadedIndices) {
  const noroBestCapsid = {};
  const noroBestPol = {};
  const rotaBestG = {};
  const rotaBestP = {};

  for (const r of localResults) {
    const parsed = parseSampleName(r.name);
    if (!parsed) continue;
    const idx = parsed.index;

    if (r.noro.capsid) {
      const score = r.noro.capsidIdentity;
      if (!noroBestCapsid[idx] || score > noroBestCapsid[idx].score) {
        noroBestCapsid[idx] = { genotype: r.noro.capsid, score };
      }
    }
    if (r.noro.polymerase) {
      const score = r.noro.polymeraseIdentity;
      if (!noroBestPol[idx] || score > noroBestPol[idx].score) {
        noroBestPol[idx] = { genotype: r.noro.polymerase, score };
      }
    }

    // For rota: V7 files → G type, V4 files → P type
    if (r.rota.vp7) {
      const score = r.rota.vp7Identity;
      if (parsed.type === 'V7' || (parsed.type !== 'V4')) {
        if (!rotaBestG[idx] || score > rotaBestG[idx].score) {
          rotaBestG[idx] = { genotype: r.rota.vp7, score };
        }
      }
    }
    if (r.rota.vp4) {
      const score = r.rota.vp4Identity;
      if (parsed.type === 'V4' || (parsed.type !== 'V7')) {
        if (!rotaBestP[idx] || score > rotaBestP[idx].score) {
          rotaBestP[idx] = { genotype: r.rota.vp4, score };
        }
      }
    }
  }

  return buildMasterTable(uploadedIndices, {
    noroBestCapsid, noroBestPol, rotaBestG, rotaBestP,
  });
}

function buildMasterTable(uploadedIndices, { noroBestCapsid, noroBestPol, rotaBestG, rotaBestP }) {
  const rotaSummary = [];
  const noroSummary = [];
  const master = [];

  const sortedIndices = [...uploadedIndices.keys()].sort((a, b) => a - b);

  // Build rota summary
  for (const idx of sortedIndices) {
    const types = uploadedIndices.get(idx);
    const g = rotaBestG[idx]?.genotype || '';
    const p = rotaBestP[idx]?.genotype || '';
    let res = '';
    if (g && p) res = `${g}${p}`;
    else if (g) res = g;
    else if (p) res = p;
    else if (types.has('V4') || types.has('V7')) res = 'failed typing';
    if (res) rotaSummary.push({ index: idx, typing_result: res });
  }

  // Build noro summary
  for (const idx of sortedIndices) {
    const types = uploadedIndices.get(idx);
    const capsid = noroBestCapsid[idx]?.genotype || '';
    const pol = noroBestPol[idx]?.genotype || '';
    let res = '';
    if (capsid || pol) {
      const pMatch = pol.match(/P\d+$/);
      const pDesignation = pMatch ? pMatch[0] : (pol || 'Could not assign');
      const capDesignation = capsid || 'Could not assign';
      res = `${capDesignation}[${pDesignation}]`;
    } else if (types.has('GI') || types.has('GII')) {
      res = 'failed typing';
    }
    if (res) noroSummary.push({ index: idx, typing_result: res });
  }

  // Build master 1-72 table
  const rotaByIndex = {};
  for (const r of rotaSummary) rotaByIndex[r.index] = r.typing_result;
  const noroByIndex = {};
  for (const r of noroSummary) noroByIndex[r.index] = r.typing_result;

  for (let i = 1; i <= 72; i++) {
    const isUploaded = uploadedIndices.has(i);
    let finalResult = '';
    if (isUploaded) {
      const parts = [];
      const rotaRes = rotaByIndex[i];
      const noroRes = noroByIndex[i];
      if (rotaRes && rotaRes !== 'failed typing') parts.push(rotaRes);
      if (noroRes && noroRes !== 'failed typing') parts.push(noroRes);
      finalResult = parts.length > 0 ? parts.join(', ') : (isUploaded ? 'failed typing' : '');
    }
    master.push({ index: i, typing_result: finalResult });
  }

  return { rotaSummary, noroSummary, master };
}

module.exports = { parseSampleName, getUploadedIndices, summarizeOnline, summarizeLocal };
