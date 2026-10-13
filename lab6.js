const now = (typeof performance !== "undefined" && performance.now)
  ? () => performance.now() : () => Date.now();
const median = (a) => { const s = a.slice().sort((x, y) => x - y), m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
function timeMedian(fn, repeats = 4) {
  const t = [];
  for (let r = 0; r < repeats; r++) { const t0 = now(); fn(); if (r > 0) t.push((now() - t0) / 1000); }
  return median(t);
}
const MB = (bytes) => (bytes / 1048576).toFixed(2) + " МБ";

function buildGrid(m, contrast) {
  const n = m * m;
  const kOf = (r, c) => (c < m / 2 ? 1 : contrast);
  const harm = (a, b) => 2 * a * b / (a + b);
  const rows = [], cols = [], vals = [];
  const diag = new Float64Array(n);
  for (let r = 0; r < m; r++) for (let c = 0; c < m; c++) {
    const i = r * m + c, k = kOf(r, c);
    for (const [dr, dc] of [[-1, 0], [1, 0], [0, -1], [0, 1]]) {
      const rr = r + dr, cc = c + dc;
      if (rr < 0 || rr >= m || cc < 0 || cc >= m) { diag[i] += k; continue; }
      const w = harm(k, kOf(rr, cc));
      diag[i] += w;
      rows.push(i); cols.push(rr * m + cc); vals.push(-w);
    }
  }
  for (let i = 0; i < n; i++) { rows.push(i); cols.push(i); vals.push(diag[i]); }
  const b = new Float64Array(n).fill(1);
  return { n, m, rows, cols, vals, b };
}

function toCSR(n, rows, cols, vals) {
  const rowPtr = new Int32Array(n + 1);
  for (const r of rows) rowPtr[r + 1]++;
  for (let i = 0; i < n; i++) rowPtr[i + 1] += rowPtr[i];
  const nnz = rows.length, colIdx = new Int32Array(nnz), val = new Float64Array(nnz);
  const next = rowPtr.slice(0, n);
  for (let t = 0; t < nnz; t++) { const p = next[rows[t]]++; colIdx[p] = cols[t]; val[p] = vals[t]; }
  return { n, nnz, rowPtr, colIdx, val };
}
function csrMul(A, x, y) {
  for (let i = 0; i < A.n; i++) {
    let s = 0;
    for (let p = A.rowPtr[i]; p < A.rowPtr[i + 1]; p++) s += A.val[p] * x[A.colIdx[p]];
    y[i] = s;
  }
}
function toDense(A) {
  const D = new Float64Array(A.n * A.n);
  for (let i = 0; i < A.n; i++) for (let p = A.rowPtr[i]; p < A.rowPtr[i + 1]; p++) D[i * A.n + A.colIdx[p]] = A.val[p];
  return D;
}
function denseMul(D, n, x, y) {
  for (let i = 0; i < n; i++) { let s = 0; const o = i * n; for (let j = 0; j < n; j++) s += D[o + j] * x[j]; y[i] = s; }
}
const norm = (v) => Math.sqrt(v.reduce((s, x) => s + x * x, 0));
function relResidual(A, x, b) {
  const y = new Float64Array(A.n); csrMul(A, x, y);
  for (let i = 0; i < A.n; i++) y[i] -= b[i];
  return norm(y) / norm(b);
}

function solveLU(Ain, n, b) {
  const A = Float64Array.from(Ain), x = Float64Array.from(b);
  for (let k = 0; k < n; k++) {
    let p = k, mx = Math.abs(A[k * n + k]);
    for (let i = k + 1; i < n; i++) { const v = Math.abs(A[i * n + k]); if (v > mx) { mx = v; p = i; } }
    if (mx === 0) throw new Error("Матриця вироджена");
    if (p !== k) {
      for (let j = 0; j < n; j++) { const t = A[k * n + j]; A[k * n + j] = A[p * n + j]; A[p * n + j] = t; }
      const t = x[k]; x[k] = x[p]; x[p] = t;
    }
    const d = A[k * n + k];
    for (let i = k + 1; i < n; i++) {
      const f = A[i * n + k] / d;
      A[i * n + k] = f;
      for (let j = k + 1; j < n; j++) A[i * n + j] -= f * A[k * n + j];
      x[i] -= f * x[k];
    }
  }
  for (let i = n - 1; i >= 0; i--) {
    let s = x[i];
    for (let j = i + 1; j < n; j++) s -= A[i * n + j] * x[j];
    x[i] = s / A[i * n + i];
  }
  return x;
}

function choleskyDense(Ain, n) {
  const L = new Float64Array(n * n);
  for (let j = 0; j < n; j++) {
    let d = Ain[j * n + j];
    for (let k = 0; k < j; k++) d -= L[j * n + k] * L[j * n + k];
    if (d <= 0) throw new Error("Матриця не додатно визначена");
    const ljj = Math.sqrt(d); L[j * n + j] = ljj;
    for (let i = j + 1; i < n; i++) {
      let s = Ain[i * n + j];
      for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
      L[i * n + j] = s / ljj;
    }
  }
  return L;
}
function solveCholeskyDense(Ain, n, b) {
  const L = choleskyDense(Ain, n), y = Float64Array.from(b);
  for (let i = 0; i < n; i++) { let s = y[i]; for (let k = 0; k < i; k++) s -= L[i * n + k] * y[k]; y[i] = s / L[i * n + i]; }
  for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k * n + i] * y[k]; y[i] = s / L[i * n + i]; }
  return y;
}

function solveBandCholesky(A, w, b) {
  const n = A.n, W = w + 1, L = new Float64Array(n * W);
  for (let i = 0; i < n; i++) for (let p = A.rowPtr[i]; p < A.rowPtr[i + 1]; p++) {
    const j = A.colIdx[p]; if (j <= i && i - j <= w) L[i * W + (i - j)] = A.val[p];
  }
  for (let i = 0; i < n; i++) {
    const j0 = Math.max(0, i - w);
    for (let j = j0; j <= i; j++) {
      let s = L[i * W + (i - j)];
      const k0 = Math.max(j0, j - w);
      for (let k = k0; k < j; k++) s -= L[i * W + (i - k)] * L[j * W + (j - k)];
      if (j === i) { if (s <= 0) throw new Error("не SPD"); L[i * W] = Math.sqrt(s); }
      else L[i * W + (i - j)] = s / L[j * W];
    }
  }
  const y = Float64Array.from(b);
  for (let i = 0; i < n; i++) { let s = y[i]; for (let k = Math.max(0, i - w); k < i; k++) s -= L[i * W + (i - k)] * y[k]; y[i] = s / L[i * W]; }
  for (let i = n - 1; i >= 0; i--) {
    let s = y[i];
    for (let k = i + 1; k <= Math.min(n - 1, i + w); k++) s -= L[k * W + (k - i)] * y[k];
    y[i] = s / L[i * W];
  }
  return y;
}

function adjacencyOf(A) {
  const adj = Array.from({ length: A.n }, () => new Set());
  for (let i = 0; i < A.n; i++) for (let p = A.rowPtr[i]; p < A.rowPtr[i + 1]; p++) {
    const j = A.colIdx[p]; if (j !== i) adj[i].add(j);
  }
  return adj;
}
function symbolicNnzL(adjIn, minDegree) {
  const n = adjIn.length, adj = adjIn.map(s => new Set(s)), done = new Uint8Array(n), perm = [];
  let nnz = 0;
  for (let step = 0; step < n; step++) {
    let v = step;
    if (minDegree) { let best = Infinity; for (let u = 0; u < n; u++) if (!done[u] && adj[u].size < best) { best = adj[u].size; v = u; } }
    done[v] = 1; perm.push(v);
    const nb = Array.from(adj[v]);
    nnz += nb.length + 1;
    for (const a of nb) {
      adj[a].delete(v);
      for (const c of nb) if (c !== a) adj[a].add(c);
    }
    adj[v].clear();
  }
  return { nnz, perm };
}
function permuteDense(D, n, perm) {
  const P = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) P[i * n + j] = D[perm[i] * n + perm[j]];
  return P;
}
function countNnzLower(L, n) { let c = 0; for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) if (Math.abs(L[i * n + j]) > 1e-14) c++; return c; }

console.log("=== 1. Коректність: нев'язка r = ||Ax-b||/||b||, m = 20, n = 400 ===");
{
  const m = 20;
  console.log("контраст | LU (півотинг)  | Холецький (щільний) | Холецький (стрічковий)");
  for (const c of [1, 10, 100, 1000]) {
    const g = buildGrid(m, c), A = toCSR(g.n, g.rows, g.cols, g.vals), D = toDense(A);
    const r1 = relResidual(A, solveLU(D, g.n, g.b), g.b);
    const r2 = relResidual(A, solveCholeskyDense(D, g.n, g.b), g.b);
    const r3 = relResidual(A, solveBandCholesky(A, m, g.b), g.b);
    console.log(`${String(c).padEnd(9)}| ${r1.toExponential(2).padEnd(15)}| ${r2.toExponential(2).padEnd(20)}| ${r3.toExponential(2)}`);
  }
}

console.log("\n=== 2. Час — пам'ять — нев'язка (контраст = 100) ===");
console.log("m   n      | LU: с / пам'ять / r           | Холецький щільн.: с / пам'ять / r | Стрічк. Холецький: с / пам'ять / r");
for (const m of [10, 15, 20, 25, 30, 50, 100, 150]) {
  const g = buildGrid(m, 100), A = toCSR(g.n, g.rows, g.cols, g.vals), n = g.n;
  let s1 = "—", s2 = "—";
  if (m <= 30) {
    const D = toDense(A);
    const tl = timeMedian(() => solveLU(D, n, g.b)), rl = relResidual(A, solveLU(D, n, g.b), g.b);
    const tc = timeMedian(() => solveCholeskyDense(D, n, g.b)), rc = relResidual(A, solveCholeskyDense(D, n, g.b), g.b);
    s1 = `${tl.toFixed(4)} / ${MB(8 * n * n)} / ${rl.toExponential(1)}`;
    s2 = `${tc.toFixed(4)} / ${MB(8 * n * n)} / ${rc.toExponential(1)}`;
  }
  const tb = timeMedian(() => solveBandCholesky(A, m, g.b), 3), rb = relResidual(A, solveBandCholesky(A, m, g.b), g.b);
  console.log(`${String(m).padEnd(4)}${String(n).padEnd(7)}| ${s1.padEnd(29)}| ${s2.padEnd(34)}| ${tb.toFixed(4)} / ${MB(8 * n * (m + 1))} / ${rb.toExponential(1)}`);
}

console.log("\n=== 3. Щільне подання проти CSR (контраст = 100) ===");
console.log("m   n      nnz     | щільне, МБ | CSR, МБ | економія | побудова CSR, мс | A·x щільне, мс | A·x CSR, мс");
for (const m of [10, 20, 30, 50, 100, 150]) {
  const g = buildGrid(m, 100), n = g.n;
  const tBuild = timeMedian(() => toCSR(n, g.rows, g.cols, g.vals), 5) * 1000;
  const A = toCSR(n, g.rows, g.cols, g.vals);
  const denseBytes = 8 * n * n, csrBytes = 8 * A.nnz + 4 * A.nnz + 4 * (n + 1);
  const x = new Float64Array(n).fill(1), y = new Float64Array(n);
  const reps = 50;
  const tCsr = timeMedian(() => { for (let r = 0; r < reps; r++) csrMul(A, x, y); }, 5) * 1000 / reps;
  let tDense = "—";
  if (m <= 50) { const D = toDense(A); tDense = (timeMedian(() => { for (let r = 0; r < 5; r++) denseMul(D, n, x, y); }, 5) * 1000 / 5).toFixed(3); }
  console.log(`${String(m).padEnd(4)}${String(n).padEnd(7)}${String(A.nnz).padEnd(8)}| ${(denseBytes / 1048576).toFixed(2).padEnd(11)}| ` +
    `${(csrBytes / 1048576).toFixed(3).padEnd(8)}| ×${(denseBytes / csrBytes).toFixed(0).padEnd(8)}| ${tBuild.toFixed(3).padEnd(17)}| ${String(tDense).padEnd(15)}| ${tCsr.toFixed(3)}`);
}

console.log("\n=== 4. Ефект впорядкування: заповнення L ===");
{
  const g = buildGrid(12, 1), A = toCSR(g.n, g.rows, g.cols, g.vals), D = toDense(A), adj = adjacencyOf(A);
  const nat = symbolicNnzL(adj, false), md = symbolicNnzL(adj, true);
  const numNat = countNnzLower(choleskyDense(D, g.n), g.n);
  const numMd = countNnzLower(choleskyDense(permuteDense(D, g.n, md.perm), g.n), g.n);
  console.log(`Перевірка (m=12): природний — символьно ${nat.nnz}, чисельно ${numNat}; мін. степінь — символьно ${md.nnz}, чисельно ${numMd}`);
}
console.log("m   n      nnz(A,нижн.) | nnz(L) природн. | nnz(L) мін.степ. | заповнення природн. / мін.степ. | у скільки разів менше");
for (const m of [20, 30, 40, 60, 80]) {
  const g = buildGrid(m, 1), A = toCSR(g.n, g.rows, g.cols, g.vals), n = g.n, adj = adjacencyOf(A);
  let lowerA = 0; for (let i = 0; i < n; i++) for (let p = A.rowPtr[i]; p < A.rowPtr[i + 1]; p++) if (A.colIdx[p] <= i) lowerA++;
  const nat = symbolicNnzL(adj, false), md = symbolicNnzL(adj, true);
  console.log(`${String(m).padEnd(4)}${String(n).padEnd(7)}${String(lowerA).padEnd(13)}| ${String(nat.nnz).padEnd(16)}| ${String(md.nnz).padEnd(17)}| ` +
    `${(nat.nnz - lowerA) + " / " + (md.nnz - lowerA)} | ×${(nat.nnz / md.nnz).toFixed(2)}`);
}
