const now = (typeof performance !== "undefined" && performance.now)
  ? () => performance.now() : () => Date.now();

function median(a) {
  const s = a.slice().sort((x, y) => x - y), m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}
function quantile(sorted, q) {
  const p = (sorted.length - 1) * q, lo = Math.floor(p), hi = Math.ceil(p);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (p - lo);
}
function bench(action, repeats = 5) {
  const t = [];
  for (let r = 0; r < repeats; r++) {
    const t0 = now(); action(); const dt = (now() - t0) / 1000;
    if (r > 0) t.push(dt);
  }
  const s = t.slice().sort((a, b) => a - b);
  return { med: median(t), min: s[0], max: s[s.length - 1], iqr: quantile(s, 0.75) - quantile(s, 0.25) };
}

function peakKernel(iters) {
  let a0 = 1, a1 = 1.1, a2 = 1.2, a3 = 1.3, a4 = 1.4, a5 = 1.5, a6 = 1.6, a7 = 1.7;
  const x = 0.9999999, y = 1e-9;
  for (let i = 0; i < iters; i++) {
    a0 = a0 * x + y; a1 = a1 * x + y; a2 = a2 * x + y; a3 = a3 * x + y;
    a4 = a4 * x + y; a5 = a5 * x + y; a6 = a6 * x + y; a7 = a7 * x + y;
  }
  return a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7;
}
const PEAK_ITERS = 20000000;
let sink = 0;
const tPeak = bench(() => { sink += peakKernel(PEAK_ITERS); }, 6);
const Ppeak = PEAK_ITERS * 8 * 2 / tPeak.med / 1e9;

const MEM_N = 4000000;
const ma = new Float64Array(MEM_N).fill(1), mb = new Float64Array(MEM_N).fill(2), mc = new Float64Array(MEM_N);
function triad(s) { for (let i = 0; i < MEM_N; i++) mc[i] = ma[i] + s * mb[i]; }
const tMem = bench(() => triad(3.0), 8);
const Bmem = (MEM_N * 3 * 8) / tMem.med / 1e9;

function heatSteps(n, steps) {
  let u = new Float64Array(n * n), v = new Float64Array(n * n);
  u[(n >> 1) * n + (n >> 1)] = 100;
  const r = 0.2;
  for (let s = 0; s < steps; s++) {
    for (let i = 1; i < n - 1; i++) {
      const row = i * n;
      for (let j = 1; j < n - 1; j++) {
        const k = row + j, uc = u[k];
        v[k] = uc + r * (u[k - n] + u[k + n] + u[k - 1] + u[k + 1] - 4 * uc);
      }
    }
    const t = u; u = v; v = t;
  }
  return u;
}

{
  const n = 41, u = heatSteps(n, 20);
  let sum = 0; for (let i = 0; i < u.length; i++) sum += u[i];
  console.log(`Перевірка: сума тепла після 20 кроків = ${sum.toFixed(6)} (очікується 100)`);
}

const FLOP_PER_CELL = 7, BYTES_PER_CELL = 16;
const I = FLOP_PER_CELL / BYTES_PER_CELL;
const Istar = Ppeak / Bmem;
const roof = Math.min(Ppeak, Bmem * I);

console.log("\n=== Еталонні константи цієї машини ===");
console.log(`P_peak = ${Ppeak.toFixed(2)} ГФлоп/с  (медіана ${tPeak.med.toFixed(4)} с, min–max ${tPeak.min.toFixed(4)}–${tPeak.max.toFixed(4)})`);
console.log(`B_mem  = ${Bmem.toFixed(2)} ГБ/с     (медіана ${tMem.med.toFixed(4)} с, min–max ${tMem.min.toFixed(4)}–${tMem.max.toFixed(4)})`);
console.log(`Точка перетину I* = P_peak/B_mem = ${Istar.toFixed(3)} FLOP/байт`);
console.log(`\nІнтенсивність ядра I = ${FLOP_PER_CELL}/${BYTES_PER_CELL} = ${I.toFixed(4)} FLOP/байт`);
console.log(`Roofline-межа для I: min(P_peak, B_mem·I) = ${roof.toFixed(2)} ГФлоп/с -> ` +
  (I < Istar ? "зона, обмежена ПАМ'ЯТТЮ" : "зона, обмежена ОБЧИСЛЕННЯМИ"));

const QUICK = false;
const STEPS = QUICK ? 100 : 500, sizes = QUICK ? [100, 150, 200, 250, 300] : [200, 300, 400, 500, 600, 800];
console.log(`\n=== Стенсиль, ${STEPS} кроків, повтори 5 (перший відкинуто) ===`);
console.log("n".padEnd(6) + "медіана,с".padEnd(12) + "min–max,с".padEnd(20) + "IQR,с".padEnd(10) + "ГФлоп/с".padEnd(10) + "% межі");
for (const n of sizes) {
  const t = bench(() => { sink += heatSteps(n, STEPS)[n + 1]; }, 5);
  const cells = (n - 2) * (n - 2) * STEPS;
  const gf = cells * FLOP_PER_CELL / t.med / 1e9;
  console.log(String(n).padEnd(6) + t.med.toFixed(4).padEnd(12) +
    `${t.min.toFixed(4)}–${t.max.toFixed(4)}`.padEnd(20) + t.iqr.toFixed(4).padEnd(10) +
    gf.toFixed(3).padEnd(10) + (100 * gf / roof).toFixed(1) + " %");
}
if (sink === 42) console.log("");
