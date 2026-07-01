/**
 * Local virus typing engine — reference-based alignment.
 * Supports Norovirus (VP1/RdRp) and Rotavirus A (VP7/VP4).
 */

const fs = require('fs');
const path = require('path');

const IUPAC = {
  A:'A', C:'C', G:'G', T:'T', U:'T',
  R:'AG', Y:'CT', S:'GC', W:'AT', K:'GT', M:'AC',
  B:'CGT', D:'AGT', H:'ACT', V:'ACG', N:'ACGT',
};
const COMPLEMENT = {
  A:'T', C:'G', G:'C', T:'A', U:'A',
  R:'Y', Y:'R', S:'S', W:'W', K:'M', M:'K',
  B:'V', D:'H', H:'D', V:'B', N:'N',
};

function basesCompatible(a, b) {
  const aSet = IUPAC[a] || '', bSet = IUPAC[b] || '';
  if (!aSet || !bSet) return false;
  return [...aSet].some(base => bSet.includes(base));
}

function reverseComplement(seq) {
  return [...seq].reverse().map(b => COMPLEMENT[b] || 'N').join('');
}

function cleanSequence(val) {
  return val.toUpperCase().replace(/[^ACGTURYSWKMBDHVN-]/g, '').replace(/U/g, 'T');
}

function parseFasta(text) {
  const records = [];
  let current = null;
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith('>')) {
      if (current) records.push(current);
      current = { name: line.slice(1).trim(), sequence: '' };
    } else if (current) {
      current.sequence += cleanSequence(line);
    }
  }
  if (current) records.push(current);
  return records.filter(r => r.sequence.length > 0);
}

function bestUngappedOverlap(query, refSeq) {
  let best = { identity: 0, matches: 0, compared: 0, overlap: 0, offset: 0 };
  const qlen = query.length, rlen = refSeq.length;
  const minUseful = Math.min(60, qlen, rlen);
  for (let offset = -rlen + minUseful; offset <= qlen - minUseful; offset++) {
    const qStart = Math.max(0, offset);
    const rStart = Math.max(0, -offset);
    const overlap = Math.min(qlen - qStart, rlen - rStart);
    if (overlap < minUseful) continue;
    let matches = 0, compared = 0;
    for (let i = 0; i < overlap; i++) {
      const qBase = query[qStart + i], rBase = refSeq[rStart + i];
      if (qBase === '-' || rBase === '-') continue;
      compared++;
      if (basesCompatible(qBase, rBase)) matches++;
    }
    const identity = compared ? (matches / compared) * 100 : 0;
    if (identity > best.identity ||
        (identity === best.identity && compared > best.compared) ||
        (identity === best.identity && compared === best.compared && overlap > best.overlap)) {
      best = { identity, matches, compared, overlap, offset };
    }
  }
  return best;
}

function bestOrientation(query, refSeq) {
  const fwd = bestUngappedOverlap(query, refSeq);
  const rev = bestUngappedOverlap(reverseComplement(query), refSeq);
  return rev.identity > fwd.identity ? { ...rev, direction: 'reverse' } : { ...fwd, direction: 'forward' };
}

function buildKmerIndex(refs, k) {
  const idx = new Map();
  for (let i = 0; i < refs.length; i++) {
    const seq = refs[i].sequence;
    const seen = new Set();
    for (let j = 0; j <= seq.length - k; j++) {
      const mer = seq.slice(j, j + k);
      if (!/^[ACGT]+$/.test(mer)) continue;
      if (seen.has(mer)) continue;
      seen.add(mer);
      if (!idx.has(mer)) idx.set(mer, []);
      idx.get(mer).push(i);
    }
  }
  return idx;
}

function rankByKmer(query, segment, refs, kmerIdx, options) {
  const k = options.k || 13;
  const candidateLimit = options.candidateLimit || 40;
  const minOverlap = options.minOverlap || 100;

  const hits = new Map();
  for (let j = 0; j <= query.length - k; j++) {
    const mer = query.slice(j, j + k);
    const entries = kmerIdx.get(mer);
    if (!entries) continue;
    for (const refIdx of entries) {
      if (refs[refIdx].segment !== segment) continue;
      hits.set(refIdx, (hits.get(refIdx) || 0) + 1);
    }
  }
  const candidates = [...hits.entries()]
    .sort((a, b) => b[1] - a[1])
    .slice(0, candidateLimit * 2)
    .map(([idx]) => idx);

  if (candidates.length === 0) {
    for (let i = 0; i < refs.length; i++) {
      if (refs[i].segment === segment) candidates.push(i);
    }
  }

  return candidates
    .map(idx => {
      const ref = refs[idx];
      const match = bestOrientation(query, ref.sequence);
      return { ...match, ref };
    })
    .filter(m => m.overlap >= Math.min(minOverlap, query.length, m.ref.sequence.length))
    .sort((a, b) => b.identity - a.identity || b.matches - a.matches || b.overlap - a.overlap)
    .slice(0, 5);
}

const NORO_REF_FILES = [
  { segment: 'RdRp', genogroup: 'GI',  file: 'reference/noro/nov-rdrp-GI.fa' },
  { segment: 'RdRp', genogroup: 'GII', file: 'reference/noro/nov-rdrp-GII.fa' },
  { segment: 'RdRp', genogroup: 'GIV', file: 'reference/noro/nov-rdrp-GIV.fa' },
  { segment: 'VP1',  genogroup: 'GI',    file: 'reference/noro/nov-vp1-GI.fa' },
  { segment: 'VP1',  genogroup: 'GII',   file: 'reference/noro/nov-vp1-GII.fa' },
  { segment: 'VP1',  genogroup: 'GIV',   file: 'reference/noro/nov-vp1-GIV.fa' },
  { segment: 'VP1',  genogroup: 'GVIII', file: 'reference/noro/nov-vp1-GVIII.fa' },
  { segment: 'VP1',  genogroup: 'GIX',   file: 'reference/noro/nov-vp1-GIX.fa' },
];

const ROTA_REF_FILES = [
  { segment: 'VP7', genotype: 'G1',  file: 'reference/rota/rva-vp7-G1.fa' },
  { segment: 'VP7', genotype: 'G2',  file: 'reference/rota/rva-vp7-G2.fa' },
  { segment: 'VP7', genotype: 'G3',  file: 'reference/rota/rva-vp7-G3.fa' },
  { segment: 'VP7', genotype: 'G4',  file: 'reference/rota/rva-vp7-G4.fa' },
  { segment: 'VP7', genotype: 'G8',  file: 'reference/rota/rva-vp7-G8.fa' },
  { segment: 'VP7', genotype: 'G9',  file: 'reference/rota/rva-vp7-G9.fa' },
  { segment: 'VP7', genotype: 'G10', file: 'reference/rota/rva-vp7-G10.fa' },
  { segment: 'VP7', genotype: 'G11', file: 'reference/rota/rva-vp7-G11.fa' },
  { segment: 'VP7', genotype: 'G12', file: 'reference/rota/rva-vp7-G12.fa' },
  { segment: 'VP4', genotype: 'P[4]',  file: 'reference/rota/rva-vp4-P[4].fa' },
  { segment: 'VP4', genotype: 'P[6]',  file: 'reference/rota/rva-vp4-P[6].fa' },
  { segment: 'VP4', genotype: 'P[8]',  file: 'reference/rota/rva-vp4-P[8].fa' },
  { segment: 'VP4', genotype: 'P[9]',  file: 'reference/rota/rva-vp4-P[9].fa' },
  { segment: 'VP4', genotype: 'P[11]', file: 'reference/rota/rva-vp4-P[11].fa' },
  { segment: 'VP4', genotype: 'P[14]', file: 'reference/rota/rva-vp4-P[14].fa' },
  { segment: 'VP4', genotype: 'P[19]', file: 'reference/rota/rva-vp4-P[19].fa' },
];

function parseNoroGenotype(name, segment) {
  const first = name.split(/\s+/)[0] || '';
  if (segment === 'RdRp') return first.match(/^G[IVX]+\.P[\w.\-]+/i)?.[0] || first;
  return first.match(/^G[IVX]+\.[\w.\-]+/i)?.[0] || first;
}

function loadReferences(baseDir) {
  const refs = [];

  for (const meta of NORO_REF_FILES) {
    const filePath = path.join(baseDir, meta.file);
    if (!fs.existsSync(filePath)) continue;
    const text = fs.readFileSync(filePath, 'utf8');
    const records = parseFasta(text);
    for (const r of records) {
      refs.push({
        ...r,
        segment: meta.segment,
        genotype: parseNoroGenotype(r.name, meta.segment),
        virus: 'noro',
        genogroup: meta.genogroup,
      });
    }
  }

  for (const meta of ROTA_REF_FILES) {
    const filePath = path.join(baseDir, meta.file);
    if (!fs.existsSync(filePath)) continue;
    const text = fs.readFileSync(filePath, 'utf8');
    const records = parseFasta(text);
    for (const r of records) {
      refs.push({
        ...r,
        segment: meta.segment,
        genotype: meta.genotype,
        virus: 'rota',
      });
    }
  }

  return refs;
}

function typeSequence(name, seq, noroRefs, rotaRefs, noroIdx, rotaIdx, options) {
  const opts = {
    k: options.k || 13,
    candidateLimit: options.candidateLimit || 40,
    minOverlap: options.minOverlap || 100,
    noroThresholds: options.noroThresholds || { VP1: 80, RdRp: 85 },
    rotaThresholds: options.rotaThresholds || { VP7: 80, VP4: 85 },
  };

  const noroVp1 = noroIdx
    ? rankByKmer(seq, 'VP1', noroRefs, noroIdx, { ...opts })
    : [];
  const noroRdrp = noroIdx
    ? rankByKmer(seq, 'RdRp', noroRefs, noroIdx, { ...opts })
    : [];
  const rotaVp7 = rotaIdx
    ? rankByKmer(seq, 'VP7', rotaRefs, rotaIdx, { ...opts })
    : [];
  const rotaVp4 = rotaIdx
    ? rankByKmer(seq, 'VP4', rotaRefs, rotaIdx, { ...opts })
    : [];

  const topNoroVp1 = noroVp1[0] || null;
  const topNoroRdrp = noroRdrp[0] || null;
  const topRotaVp7 = rotaVp7[0] || null;
  const topRotaVp4 = rotaVp4[0] || null;

  const noroVp1Pass = topNoroVp1 && topNoroVp1.identity >= opts.noroThresholds.VP1 && topNoroVp1.overlap >= opts.minOverlap;
  const noroRdrpPass = topNoroRdrp && topNoroRdrp.identity >= opts.noroThresholds.RdRp && topNoroRdrp.overlap >= opts.minOverlap;
  const rotaVp7Pass = topRotaVp7 && topRotaVp7.identity >= opts.rotaThresholds.VP7 && topRotaVp7.overlap >= opts.minOverlap;
  const rotaVp4Pass = topRotaVp4 && topRotaVp4.identity >= opts.rotaThresholds.VP4 && topRotaVp4.overlap >= opts.minOverlap;

  const noroBestId = Math.max(topNoroVp1?.identity || 0, topNoroRdrp?.identity || 0);
  const rotaBestId = Math.max(topRotaVp7?.identity || 0, topRotaVp4?.identity || 0);
  const noroPassCount = (noroVp1Pass ? 1 : 0) + (noroRdrpPass ? 1 : 0);
  const rotaPassCount = (rotaVp7Pass ? 1 : 0) + (rotaVp4Pass ? 1 : 0);

  return {
    name,
    length: seq.length,
    noro: {
      capsid: noroVp1Pass ? topNoroVp1.ref.genotype : '',
      capsidIdentity: topNoroVp1?.identity || 0,
      capsidOverlap: topNoroVp1?.overlap || 0,
      polymerase: noroRdrpPass ? topNoroRdrp.ref.genotype : '',
      polymeraseIdentity: topNoroRdrp?.identity || 0,
      polymeraseOverlap: topNoroRdrp?.overlap || 0,
    },
    rota: {
      vp7: rotaVp7Pass ? topRotaVp7.ref.genotype : '',
      vp7Identity: topRotaVp7?.identity || 0,
      vp7Overlap: topRotaVp7?.overlap || 0,
      vp4: rotaVp4Pass ? topRotaVp4.ref.genotype : '',
      vp4Identity: topRotaVp4?.identity || 0,
      vp4Overlap: topRotaVp4?.overlap || 0,
    },
    bestVirus: noroPassCount >= rotaPassCount && noroBestId >= rotaBestId ? 'noro' : 'rota',
    noroBestId,
    rotaBestId,
    noroPassCount,
    rotaPassCount,
  };
}

function generateNorovirusCSV(results) {
  const lines = ['"Name","Length","Family Genus Genogroup","BLAST score","polymerase_1","polymerase_2","capsid_1","capsid_2","Report","genome"'];
  for (const r of results) {
    const pol = r.noro.polymerase || '';
    const capsid = r.noro.capsid || '';
    const identity = Math.max(r.noro.capsidIdentity, r.noro.polymeraseIdentity).toFixed(2);
    lines.push(`"${r.name}",${r.length},"Caliciviridae Norovirus",${identity},"${pol}","","${capsid}","","Report",""`);
  }
  return lines.join('\n');
}

function generateRotavirusCSV(results) {
  const lines = ['"Name","Length","BLAST result","BLAST score","cluster/support_1","cluster/support_2","Report","genome"'];
  for (const r of results) {
    const hasG = !!r.rota.vp7;
    const hasP = !!r.rota.vp4;
    let cluster = '', support = '';
    if (hasG && hasP) { cluster = r.rota.vp7; support = r.rota.vp4; }
    else if (hasG) { cluster = r.rota.vp7; }
    else if (hasP) { cluster = r.rota.vp4; }
    const bestId = Math.max(r.rota.vp7Identity, r.rota.vp4Identity).toFixed(2);
    lines.push(`"${r.name}",${r.length},"Reoviridae Rotavirus A",${bestId},"${cluster}","${support}","Report",""`);
  }
  return lines.join('\n');
}

module.exports = {
  parseFasta, cleanSequence, loadReferences, buildKmerIndex,
  typeSequence, generateNorovirusCSV, generateRotavirusCSV,
  bestOrientation, NORO_REF_FILES, ROTA_REF_FILES,
};
