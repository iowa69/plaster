/**
 * api.js — thin, defensive client for the Assemblage HTTP API.
 *
 * Every call resolves to parsed JSON or throws an `ApiError` carrying the
 * server's `error` string verbatim (plus `detail`, which the UI shows in a
 * smaller font). Network failures are turned into ApiErrors too, so callers
 * only ever have to handle one error type.
 */

export class ApiError extends Error {
  constructor(message, detail, status, endpoint) {
    super(message || 'Request failed');
    this.name = 'ApiError';
    /** The server's `error` field, shown to the user verbatim. */
    this.error = message || 'Request failed';
    this.detail = detail || '';
    this.status = status || 0;
    this.endpoint = endpoint || '';
  }
}

/**
 * Work out where `/api` lives relative to the page we were served from, so the
 * UI keeps working if the server ever mounts the static files under a prefix.
 */
function computeBase() {
  try {
    let p = window.location.pathname || '/';
    // Strip a trailing document name (index.html or anything with a dot).
    const slash = p.lastIndexOf('/');
    const tail = p.slice(slash + 1);
    if (tail.includes('.')) p = p.slice(0, slash + 1);
    if (!p.endsWith('/')) p += '/';
    return p + 'api';
  } catch {
    return '/api';
  }
}

export const API_BASE = computeBase();

/** Build a query string, skipping null/undefined/'' values. */
function qs(params) {
  if (!params) return '';
  const sp = new URLSearchParams();
  for (const [k, v] of Object.entries(params)) {
    if (v === null || v === undefined || v === '') continue;
    sp.append(k, String(v));
  }
  const s = sp.toString();
  return s ? '?' + s : '';
}

async function request(method, path, opts = {}) {
  const url = API_BASE + path + qs(opts.query);
  const init = { method, headers: { Accept: 'application/json' } };
  if (opts.signal) init.signal = opts.signal;
  if (opts.body !== undefined) {
    init.headers['Content-Type'] = 'application/json';
    init.body = JSON.stringify(opts.body);
  }

  let res;
  try {
    res = await fetch(url, init);
  } catch (e) {
    if (e && e.name === 'AbortError') throw e;
    throw new ApiError(
      'Cannot reach the Assemblage server',
      String((e && e.message) || e) + '\n' + method + ' ' + url,
      0, path,
    );
  }

  const text = await res.text();
  let data = null;
  if (text) {
    try { data = JSON.parse(text); } catch { data = null; }
  }

  if (!res.ok) {
    const msg = (data && (data.error || data.message)) || `HTTP ${res.status} ${res.statusText || ''}`.trim();
    const detail = (data && data.detail) || (data ? '' : text.slice(0, 600));
    throw new ApiError(msg, detail, res.status, path);
  }

  // A JSON body of literal `null` is legitimate (e.g. GET /api/scaffold).
  if (data === null && text.trim() && text.trim() !== 'null') return text;
  return data;
}

const get = (p, query, signal) => request('GET', p, { query, signal });
const post = (p, body, signal) => request('POST', p, { body, signal });
const put = (p, body, signal) => request('PUT', p, { body, signal });
const del = (p, signal) => request('DELETE', p, { signal });

export const api = {
  ApiError,

  /* ---- project ---- */
  status: (signal) => get('/status', null, signal),
  load: (path, format, signal) => post('/load', { path, format: format || null }, signal),
  loadPaths: (path, signal) => post('/load-paths', { path }, signal),

  /* ---- graph ---- */
  graph: ({ min_length = 0, max_nodes = 15000, component = null } = {}, signal) =>
    get('/graph', { min_length, max_nodes, component }, signal),
  segment: (name, signal) => get('/segment/' + encodeURIComponent(name), null, signal),

  /* ---- reference ---- */
  addReference: (body, signal) => post('/reference', body, signal),
  dropReference: (signal) => del('/reference', signal),

  /* ---- metrics ---- */
  report: (signal) => get('/report', null, signal),

  /* ---- mutation ---- */
  op: (op, args, signal) => post('/op', { op, args: args || {} }, signal),
  undo: (signal) => post('/undo', {}, signal),

  /* ---- scaffolding ---- */
  buildScaffold: (body, signal) => post('/scaffold', body, signal),
  getScaffold: (signal) => get('/scaffold', null, signal),
  putScaffold: (plan, signal) => put('/scaffold', plan, signal),

  /* ---- search ---- */
  search: (query, min_identity, signal) => post('/search', { query, min_identity }, signal),

  /* ---- export ---- */
  exportFiles: (outdir, what, signal) => post('/export', { outdir, what }, signal),
  downloadUrl: (kind) => API_BASE + '/download/' + encodeURIComponent(kind),

  /* ---- session ---- */
  getSession: (signal) => get('/session', null, signal),
  saveSession: (payload, signal) => post('/session', payload, signal),

  /* ---- file picker ---- */
  browse: (path, signal) => get('/browse', { path }, signal),
};

export default api;
