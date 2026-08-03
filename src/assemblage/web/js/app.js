/**
 * app.js — entry point.
 *
 * index.html loads only this module; everything else is imported from here so
 * there is one place to look for the boot sequence.
 */

import { start } from './ui.js';

function fatal(message) {
  const box = document.getElementById('toasts') || document.body;
  const el = document.createElement('div');
  el.className = 'toast toast-error';
  el.textContent = message;
  box.appendChild(el);
  console.error(message);
}

function boot() {
  start().catch((err) => {
    fatal(`Assemblage failed to start: ${err && err.message ? err.message : err}`);
  });
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', boot, { once: true });
} else {
  boot();
}
