/**
 * Summarize typing results into a 1-72 index table.
 * Produces a combined comparison: index | local_result | online_result
 */

const fs = require('fs');
const path = require('path');

// Parse sample name: May05V4 → { month: 'May', index: 5, type: 'V4' }
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

// Scan uploaded files to discover indices and their types
function getUploadedIndices(seqDir) {
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

// Extract best genotypes per index from RIVM online results
function extractOnlineByIndex(noroResults, rotaResults) {
  const noroCapsid = {};  // idx -> {genotype, score}
  const noroPol = {};     // idx -> {genotype, score}
  const rotaG = {};       // idx -> {genotype, score}
  const rotaP = {};       // idx -> {genotype, score}

  if (noroResults) {
    for (const row of noroResults) {
      const name = row['Name'] || row['name'] || '';
      const parsed = parseSampleName(name);
      if (!parsed) continue;
      const idx = parsed.index;
      const score = parseFloat(row['BLAST score'] || 0);
      const capsid = row['capsid type'] || row['capsid_1'] || '';
      const pol = row['polymerase type'] || row['polymerase_1'] || '';
      if (capsid && (!noroCapsid[idx] || score > noroCapsid[idx].score))
        noroCapsid[idx] = { genotype: capsid, score };
      if (pol && (!noroPol[idx] || score > noroPol[idx].score))
        noroPol[idx] = { genotype: pol, score };
    }
  }

  if (rotaResults) {
    for (const row of rotaResults) {
      const name = row['Name'] || row['name'] || '';
      const parsed = parseSampleName(name);
      if (!parsed) continue;
      const idx = parsed.index;
      const score = parseFloat(row['BLAST score'] || 0);
      const val1 = row['cluster result'] || row['cluster/support_1'] || row['cluster/support type'] || '';
      const val2 = row['cluster support'] || row['cluster/support_2'] || row['cluster/support subtype'] || '';
      let gGen = '', pGen = '';
      for (const val of [val1, val2]) {
        if (!val) continue;
        if (/^G\d+$/i.test(val)) gGen = val;
        if (/^P\[\d+\]$/i.test(val)) pGen = val;
      }
      if (gGen && (!rotaG[idx] || score > rotaG[idx].score))
        rotaG[idx] = { genotype: gGen, score };
      if (pGen && (!rotaP[idx] || score > rotaP[idx].score))
        rotaP[idx] = { genotype: pGen, score };
    }
  }

  return { noroCapsid, noroPol, rotaG, rotaP };
}

// Extract best genotypes per index from local typing results
function extractLocalByIndex(localResults) {
  const noroCapsid = {};
  const noroPol = {};
  const rotaG = {};
  const rotaP = {};

  for (const r of localResults) {
    const parsed = parseSampleName(r.name);
    if (!parsed) continue;
    const idx = parsed.index;

    if (r.noro.capsid) {
      const score = r.noro.capsidIdentity;
      if (!noroCapsid[idx] || score > noroCapsid[idx].score)
        noroCapsid[idx] = { genotype: r.noro.capsid, score };
    }
    if (r.noro.polymerase) {
      const score = r.noro.polymeraseIdentity;
      if (!noroPol[idx] || score > noroPol[idx].score)
        noroPol[idx] = { genotype: r.noro.polymerase, score };
    }
    if (r.rota.vp7) {
      const score = r.rota.vp7Identity;
      if (!rotaG[idx] || score > rotaG[idx].score)
        rotaG[idx] = { genotype: r.rota.vp7, score };
    }
    if (r.rota.vp4) {
      const score = r.rota.vp4Identity;
      if (!rotaP[idx] || score > rotaP[idx].score)
        rotaP[idx] = { genotype: r.rota.vp4, score };
    }
  }

  return { noroCapsid, noroPol, rotaG, rotaP };
}

// Build a typing result string for an index from extracted genotype maps
function buildResultString(idx, types, { noroCapsid, noroPol, rotaG, rotaP }) {
  const parts = [];

  // Rotavirus: G + P
  const g = rotaG[idx]?.genotype || '';
  const p = rotaP[idx]?.genotype || '';
  if (g || p) {
    if (g && p) parts.push(`${g}${p}`);
    else if (g) parts.push(g);
    else parts.push(p);
  } else if (types && (types.has('V4') || types.has('V7'))) {
    parts.push('rota failed');
  }

  // Norovirus: capsid[P]
  const capsid = noroCapsid[idx]?.genotype || '';
  const pol = noroPol[idx]?.genotype || '';
  if (capsid || pol) {
    const pMatch = pol.match(/P[\w.\-]+$/);
    const pDesignation = pMatch ? pMatch[0] : (pol || '?');
    parts.push(`${capsid || '?'}[${pDesignation}]`);
  } else if (types && (types.has('GI') || types.has('GII'))) {
    parts.push('noro failed');
  }

  return parts.join(', ');
}

/**
 * Build combined 1-72 summary with local + online + BLAST columns.
 * Returns { master: [{index, local, online, blast}], ... }
 */
function buildComparisonSummary(uploadedIndices, localResults, onlineNoroResults, onlineRotaResults, blastResults) {
  const local = extractLocalByIndex(localResults);
  const online = extractOnlineByIndex(onlineNoroResults, onlineRotaResults);

  // Extract BLAST results by index
  const blastByIndex = {};
  if (blastResults) {
    for (const r of blastResults) {
      const parsed = parseSampleName(r.name);
      if (!parsed) continue;
      const idx = parsed.index;
      // Format: capsid[P] for noro, or G+P for rota
      const parts = [];
      if (r.rota?.gType || r.rota?.pType) {
        parts.push(`${r.rota.gType || ''}${r.rota.pType || ''}`);
      }
      if (r.noro?.capsid || r.noro?.polymerase) {
        parts.push(`${r.noro.capsid || '?'}[${r.noro.polymerase || '?'}]`);
      }
      blastByIndex[idx] = parts.join(', ') || (r.blast?.error ? `error: ${r.blast.error}` : '');
    }
  }

  const master = [];
  for (let i = 1; i <= 72; i++) {
    const types = uploadedIndices.get(i);
    const isUploaded = uploadedIndices.has(i);

    let localResult = '';
    let onlineResult = '';
    let blastResult = '';

    if (isUploaded) {
      localResult = buildResultString(i, types, local) || '';
      onlineResult = buildResultString(i, types, online) || '';
      blastResult = blastByIndex[i] || '';
    }

    master.push({ index: i, local: localResult, online: onlineResult, blast: blastResult });
  }

  return { master, local, online };
}

// Legacy wrappers for backward compatibility
function summarizeOnline(noroResults, rotaResults, uploadedIndices) {
  const data = extractOnlineByIndex(noroResults, rotaResults);
  const rotaSummary = [], noroSummary = [], master = [];
  const sorted = [...uploadedIndices.keys()].sort((a, b) => a - b);

  for (const idx of sorted) {
    const g = data.rotaG[idx]?.genotype || '';
    const p = data.rotaP[idx]?.genotype || '';
    let res = g && p ? `${g}${p}` : g || p || '';
    if (!res && (uploadedIndices.get(idx).has('V4') || uploadedIndices.get(idx).has('V7'))) res = 'failed typing';
    if (res) rotaSummary.push({ index: idx, typing_result: res });
  }

  for (const idx of sorted) {
    const capsid = data.noroCapsid[idx]?.genotype || '';
    const pol = data.noroPol[idx]?.genotype || '';
    let res = '';
    if (capsid || pol) {
      const pMatch = pol.match(/P\d+$/);
      res = `${capsid || '?'}[${pMatch ? pMatch[0] : (pol || '?')}]`;
    }
    if (res) noroSummary.push({ index: idx, typing_result: res });
  }

  for (let i = 1; i <= 72; i++) {
    const isUploaded = uploadedIndices.has(i);
    let finalResult = '';
    if (isUploaded) {
      const rotaRes = rotaSummary.find(r => r.index === i)?.typing_result;
      const noroRes = noroSummary.find(r => r.index === i)?.typing_result;
      const parts = [];
      if (rotaRes && rotaRes !== 'failed typing') parts.push(rotaRes);
      if (noroRes && noroRes !== 'failed typing') parts.push(noroRes);
      finalResult = parts.join(', ') || (isUploaded ? 'failed typing' : '');
    }
    master.push({ index: i, typing_result: finalResult });
  }

  return { rotaSummary, noroSummary, master };
}

function summarizeLocal(localResults, uploadedIndices) {
  return summarizeOnline(null, null, uploadedIndices); // stub — use buildComparisonSummary instead
}

module.exports = { parseSampleName, getUploadedIndices, buildComparisonSummary, summarizeOnline, summarizeLocal };
