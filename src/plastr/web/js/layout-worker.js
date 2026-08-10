/**
 * layout-worker.js — module worker hosting the layout engine.
 *
 * Loaded as:  new Worker('js/layout-worker.js', { type: 'module' })
 *
 * Protocol
 * --------
 * main -> worker
 *   {cmd: 'init',      data}                structural arrays (transferred)
 *   {cmd: 'positions', px, py}              overwrite coordinates
 *   {cmd: 'run',       mode, subset, params}
 *   {cmd: 'stop'}
 *
 * worker -> main
 *   {type: 'tick',  px, py, iter, alpha, mode}
 *   {type: 'done',  px, py, iter, alpha, mode}
 *   {type: 'error', message}
 *
 * The iteration loop deliberately runs in short slices separated by
 * `setTimeout(0)`: a worker is single-threaded, so a long synchronous loop
 * would never see the `stop` message and the Stop button would do nothing.
 */

import { LayoutEngine } from './layout.js';

const engine = new LayoutEngine();
let running = false;
let timer = 0;
/** Identifies the current run so stale results can be discarded. */
let runId = 0;
/** Milliseconds of computation per slice before yielding to the message queue. */
const SLICE_MS = 12;

function now() {
  return (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
}

function emit(type, res) {
  // Send copies so the engine keeps ownership of its own arrays.
  const px = Float32Array.from(engine.px);
  const py = Float32Array.from(engine.py);
  self.postMessage({
    type,
    px, py,
    iter: res.iter,
    alpha: res.alpha,
    mode: engine.mode,
    // Echoed so the controller can discard messages from a run it has
    // already replaced; without it the previous run's 'done' arrives after
    // the new one starts and reports the layout as finished.
    runId,
  }, [px.buffer, py.buffer]);
}

function loop() {
  timer = 0;
  if (!running) return;
  let res = { done: false, iter: engine.iter, alpha: engine.alpha };
  const t0 = now();
  try {
    do {
      res = engine.step();
    } while (!res.done && (now() - t0) < SLICE_MS);
  } catch (err) {
    running = false;
    self.postMessage({ type: 'error', message: 'layout failed: ' + ((err && err.message) || err) });
    return;
  }
  if (res.done) {
    running = false;
    emit('done', res);
  } else {
    emit('tick', res);
    timer = setTimeout(loop, 0);
  }
}

function stop() {
  running = false;
  engine.stop();
  if (timer) { clearTimeout(timer); timer = 0; }
}

self.onmessage = (ev) => {
  const msg = ev && ev.data;
  if (!msg || !msg.cmd) return;
  try {
    switch (msg.cmd) {
      case 'init':
        stop();
        if (!engine.init(msg.data || {})) {
          self.postMessage({ type: 'error', message: 'layout: empty graph' });
        }
        break;

      case 'positions':
        if (msg.px && msg.py) engine.setPositions(msg.px, msg.py);
        break;

      case 'run': {
        stop();
        runId = Number(msg.runId) || runId + 1;
        if (!engine.ready) {
          self.postMessage({ type: 'error', message: 'layout: no graph uploaded yet' });
          break;
        }
        const ok = engine.begin(msg.mode || 'force', msg.subset || null, msg.params || {});
        if (!ok) {
          self.postMessage({ type: 'error', message: 'layout: nothing to arrange' });
          break;
        }
        running = true;
        loop();
        break;
      }

      case 'stop':
        stop();
        if (Number.isFinite(msg.runId)) runId = Number(msg.runId);
        if (engine.ready) emit('done', { iter: engine.iter, alpha: 0 });
        else self.postMessage({ type: 'done', iter: 0, alpha: 0, mode: engine.mode, runId });
        break;

      default:
        break;
    }
  } catch (err) {
    running = false;
    self.postMessage({ type: 'error', message: 'layout worker: ' + ((err && err.message) || err) });
  }
};

// Tell the host we are alive; harmless if nobody listens.
self.postMessage({ type: 'ready' });
