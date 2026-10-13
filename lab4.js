const now = (typeof performance !== "undefined" && performance.now)
  ? () => performance.now() : () => Date.now();

function makeRng(seed) {
  let s = seed >>> 0;
  return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
}
function randomMatrix(n, rng) {
  const m = new Float64Array(n * n);
  for (let i = 0; i < m.length; i++) m[i] = rng() - 0.5;
  return m;
}

function toRows(flat, n) {
  const rows = new Array(n);
  for (let i = 0; i < n; i++) rows[i] = Array.from(flat.subarray(i * n, (i + 1) * n));
  return rows;
}
function mulNaive(A, B, n) {
  const C = new Array(n);
  for (let i = 0; i < n; i++) {
    const row = new Array(n);
    for (let j = 0; j < n; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += A[i][k] * B[k][j];
      row[j] = s;
    }
    C[i] = row;
  }
  return C;
}

function mulBlocked(A, B, n, bs = 48) {
  const C = new Float64Array(n * n);
  for (let ii = 0; ii < n; ii += bs) {
    const iMax = Math.min(ii + bs, n);
    for (let kk = 0; kk < n; kk += bs) {
      const kMax = Math.min(kk + bs, n);
      for (let jj = 0; jj < n; jj += bs) {
        const jMax = Math.min(jj + bs, n);
        for (let i = ii; i < iMax; i++) {
          for (let k = kk; k < kMax; k++) {
            const a = A[i * n + k];
            const bRow = k * n, cRow = i * n;
            for (let j = jj; j < jMax; j++) C[cRow + j] += a * B[bRow + j];
          }
        }
      }
    }
  }
  return C;
}

function mulTuned(A, B, n) {
  const Bt = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Bt[j * n + i] = B[i * n + j];
  const C = new Float64Array(n * n);
  const n2 = n - (n % 2);
  for (let i = 0; i < n2; i += 2) {
    const a0 = i * n, a1 = a0 + n;
    for (let j = 0; j < n2; j += 2) {
      const b0 = j * n, b1 = b0 + n;
      let c00 = 0, c01 = 0, c10 = 0, c11 = 0;
      for (let k = 0; k < n; k++) {
        const x0 = A[a0 + k], x1 = A[a1 + k], y0 = Bt[b0 + k], y1 = Bt[b1 + k];
        c00 += x0 * y0; c01 += x0 * y1; c10 += x1 * y0; c11 += x1 * y1;
      }
      C[i * n + j] = c00; C[i * n + j + 1] = c01;
      C[(i + 1) * n + j] = c10; C[(i + 1) * n + j + 1] = c11;
    }
  }
  for (let i = 0; i < n; i++) {
    for (let j = (i < n2 ? n2 : 0); j < n; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += A[i * n + k] * Bt[j * n + k];
      C[i * n + j] = s;
    }
  }
  return C;
}

function median(a) {
  const s = a.slice().sort((x, y) => x - y), m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}
function timeMedian(action, repeats = 5) {
  const times = [];
  for (let r = 0; r < repeats; r++) {
    const t0 = now();
    action();
    const dt = (now() - t0) / 1000;
    if (r === 0) continue;
    times.push(dt);
  }
  return median(times);
}

const T975 = [12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
  2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
  2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042];
function fitPowerLaw(xs, ys) {
  const n = xs.length;
  const lx = xs.map(Math.log), ly = ys.map(Math.log);
  const mx = lx.reduce((a, b) => a + b) / n, my = ly.reduce((a, b) => a + b) / n;
  let sxx = 0, sxy = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxx += (lx[i] - mx) ** 2; sxy += (lx[i] - mx) * (ly[i] - my); syy += (ly[i] - my) ** 2;
  }
  const b = sxy / sxx, a = my - b * mx;
  let sse = 0;
  for (let i = 0; i < n; i++) sse += (ly[i] - (a + b * lx[i])) ** 2;
  const df = n - 2;
  const se = Math.sqrt(sse / df / sxx);
  const t = df <= 30 ? T975[df - 1] : 1.96;
  return { b, lnC: a, se, lo: b - t * se, hi: b + t * se, r2: 1 - sse / syy, t, df };
}

function asciiPlot(sizes, series, names, marks) {
  const W = 60, H = 18;
  const xs = sizes.map(Math.log);
  const allY = series.flat().map(Math.log);
  const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...allY), y1 = Math.max(...allY);
  const grid = Array.from({ length: H }, () => Array(W).fill(" "));
  series.forEach((ys, si) => ys.forEach((y, i) => {
    const cx = Math.round((xs[i] - x0) / (x1 - x0) * (W - 1));
    const cy = H - 1 - Math.round((Math.log(y) - y0) / (y1 - y0) * (H - 1));
    grid[cy][cx] = marks[si];
  }));
  console.log(`ln T (c) від ${y0.toFixed(2)} до ${y1.toFixed(2)}; ln n від ${x0.toFixed(2)} до ${x1.toFixed(2)}`);
  grid.forEach(r => console.log("|" + r.join("")));
  console.log("+" + "-".repeat(W));
  console.log(names.map((n, i) => `${marks[i]} = ${n}`).join("   "));
}

const sizes = [50, 75, 100, 150, 200, 250, 300, 400];
const REPEATS = 5;
const rng = makeRng(2026);

const impl = [
  { name: "Наївна, Θ(n³)", theory: 3.0, times: [] },
  { name: "Блокова, Θ(n³)", theory: 3.0, times: [] },
  { name: "Оптимізоване ядро, Θ(n³)", theory: 3.0, times: [] },
];

console.log("Перевірка коректності (макс. |різниця| відносно наївної):");
for (const n of sizes) {
  const A = randomMatrix(n, rng), B = randomMatrix(n, rng);
  const Ar = toRows(A, n), Br = toRows(B, n);
  const ref = mulNaive(Ar, Br, n);
  const c2 = mulBlocked(A, B, n), c3 = mulTuned(A, B, n);
  let d2 = 0, d3 = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    d2 = Math.max(d2, Math.abs(ref[i][j] - c2[i * n + j]));
    d3 = Math.max(d3, Math.abs(ref[i][j] - c3[i * n + j]));
  }
  console.log(`  n=${n}: блокова ${d2.toExponential(1)}, ядро ${d3.toExponential(1)}`);
  if (d2 > 1e-9 || d3 > 1e-9) throw new Error("Реалізації дають різні результати!");

  impl[0].times.push(timeMedian(() => mulNaive(Ar, Br, n), REPEATS));
  impl[1].times.push(timeMedian(() => mulBlocked(A, B, n), REPEATS));
  impl[2].times.push(timeMedian(() => mulTuned(A, B, n), REPEATS));
}

console.log("\nЧас (медіана з " + (REPEATS - 1) + " вимірів, с):");
console.log("n".padEnd(6) + impl.map(x => x.name.padEnd(28)).join(""));
sizes.forEach((n, i) =>
  console.log(String(n).padEnd(6) + impl.map(x => x.times[i].toFixed(5).padEnd(28)).join("")));

console.log("\nОцінка показника степеня b (T = C·n^b), 95 % ДІ (t-розподіл, df = " + (sizes.length - 2) + "):");
for (const im of impl) {
  const f = fitPowerLaw(sizes, im.times);
  const inside = f.lo <= im.theory && im.theory <= f.hi;
  console.log(`${im.name}: b = ${f.b.toFixed(2)}  [${f.lo.toFixed(2)}; ${f.hi.toFixed(2)}]  ` +
    `R² = ${f.r2.toFixed(3)}  теор. = ${im.theory}  -> ` +
    (inside ? "теорія потрапляє в ДІ" : "теорія НЕ потрапляє в ДІ"));
}

console.log("\nГрафік у логарифмічних осях:");
asciiPlot(sizes, impl.map(x => x.times), impl.map(x => x.name), ["N", "B", "K"]);

console.log("\nПрискорення відносно наївної на найбільшому n=" + sizes[sizes.length - 1] + ": блокова ×" +
  (impl[0].times.at(-1) / impl[1].times.at(-1)).toFixed(1) + ", ядро ×" +
  (impl[0].times.at(-1) / impl[2].times.at(-1)).toFixed(1));
