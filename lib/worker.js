/**
 * Worker thread for local virus typing.
 * Receives { name, sequence, options } via postMessage.
 * Returns typing result via parentPort.
 */

const { parentPort, workerData } = require('worker_threads');
const {
  loadReferences, buildKmerIndex, typeSequence,
} = require('./local-typing');

// Load references once at startup
const allRefs = loadReferences(workerData.baseDir);
const noroRefs = allRefs.filter(r => r.virus === 'noro');
const rotaRefs = allRefs.filter(r => r.virus === 'rota');
const noroIdx = buildKmerIndex(noroRefs, 13);
const rotaIdx = buildKmerIndex(rotaRefs, 13);

parentPort.on('message', (task) => {
  try {
    const result = typeSequence(
      task.name, task.sequence,
      noroRefs, rotaRefs, noroIdx, rotaIdx,
      task.options,
    );
    parentPort.postMessage({ status: 'success', result });
  } catch (err) {
    parentPort.postMessage({ status: 'error', error: err.message });
  }
});
