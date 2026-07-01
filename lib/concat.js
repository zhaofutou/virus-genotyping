/**
 * Concatenate sequence files from a directory into a single FASTA.
 */

const fs = require('fs');
const path = require('path');

function concatSequences(seqDir, outputDir) {
  const allFiles = [];
  function scan(dir) {
    if (!fs.existsSync(dir)) return;
    for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
      const full = path.join(dir, entry.name);
      if (entry.isDirectory()) {
        scan(full);
      } else {
        const ext = path.extname(entry.name).toLowerCase();
        if (['.seq', '.fasta', '.fa'].includes(ext) && !entry.name.startsWith('concated')) {
          const baseName = path.basename(entry.name, ext);
          allFiles.push({ path: full, ext, baseName });
        }
      }
    }
  }
  scan(seqDir);

  if (allFiles.length === 0) {
    throw new Error('No .seq or .fasta files found in upload');
  }

  // Deduplicate: group by baseName and prefer .fasta/.fa over .seq
  const grouped = {};
  for (const f of allFiles) {
    const existing = grouped[f.baseName];
    if (!existing) {
      grouped[f.baseName] = f;
    } else {
      if (['.fasta', '.fa'].includes(f.ext) && existing.ext === '.seq') {
        grouped[f.baseName] = f;
      }
    }
  }
  const seqFiles = Object.values(grouped);

  const outPath = path.join(outputDir, 'concated2.fasta');
  let out = '';
  for (const f of seqFiles) {
    const header = f.baseName.replace(/[^a-zA-Z0-9_-]/g, '_');
    const text = fs.readFileSync(f.path, 'utf8');
    const seq = text.split('\n')
      .filter(line => !line.startsWith('>'))
      .join('')
      .replace(/\s/g, '');
    out += `>${header}\n${seq}\n`;
  }
  fs.writeFileSync(outPath, out);
  return { path: outPath, count: seqFiles.length };
}

module.exports = { concatSequences };
