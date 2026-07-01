/**
 * Concatenate sequence files from a directory into a single FASTA.
 */

const fs = require('fs');
const path = require('path');

function concatSequences(seqDir, outputDir) {
  const seqFiles = [];
  function scan(dir) {
    if (!fs.existsSync(dir)) return;
    for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
      const full = path.join(dir, entry.name);
      if (entry.isDirectory()) { scan(full); }
      else if (/\.(seq|fasta|fa)$/i.test(entry.name) && !entry.name.startsWith('concated')) {
        seqFiles.push(full);
      }
    }
  }
  scan(seqDir);

  if (seqFiles.length === 0) {
    throw new Error('No .seq or .fasta files found in upload');
  }

  const outPath = path.join(outputDir, 'concated2.fasta');
  let out = '';
  for (const f of seqFiles) {
    const fname = path.basename(f, path.extname(f));
    const header = fname.replace(/\s+/g, '_');
    const text = fs.readFileSync(f, 'utf8');
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
