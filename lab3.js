const card = {
  title: "Підбір параметрів схеми заміщення асинхронного двигуна за осцилограмою пускового струму",
  decisionVariables: "R = Rs + Rr' (Ом, неперервна, >0), L = Ls + Lr' (Гн, неперервна, >0); " +
                     "додатково (повна схема): Rs, Rr', Xs, Xr', Xm",
  criterion: "Один критерій: мінімізувати суму квадратів відхилень SSE = Σ (i_мод(t_k) − i_осц(t_k))²",
  hardConstraints: "R > 0, L > 0; ф-ла моделі виконується точно (рівняння кола)",
  softConstraints: "Параметри в фізично правдоподібних межах (штраф за вихід, напр. R ∈ [0.05; 20] Ом)",
  availableData: "Відліки струму осцилограми (N ≈ 500–5000), крок ≈ 0.1–1 мс, шум ≈ 1–3 % від амплітуди; U, f відомі",
  expectedSize: "2 параметри (спрощена модель) до ≈5 (повна схема); 10^3–10^4 відліків",
  problemType: "Ідентифікаційна (обернена) задача: відновлення параметрів моделі за вимірюваннями (табл. 1.1)",
  complexity: "Нелінійна МНК-задача з малою кількістю змінних: клас P (локальний розв'язок за поліноміальний час, " +
              "O(N) на обчислення цільової функції); NP-складність відсутня, але можливі локальні мінімуми"
};
console.log("=== КАРТКА ЗАДАЧІ ===");
for (const [k, v] of Object.entries(card)) console.log(k + ": " + v);

function current(t, R, L, Um, w) {
  const X = w * L, Z = Math.hypot(R, X), phi = Math.atan2(X, R);
  return (Um / Z) * (Math.sin(w * t - phi) + Math.sin(phi) * Math.exp(-t * R / L));
}

function makeRng(seed) {
  let s = seed >>> 0;
  return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
}
function gauss(rng) {
  const u = 1 - rng(), v = rng();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

const Um = 325, f = 50, w = 2 * Math.PI * f;
const trueR = 1.2, trueL = 0.012;
const N = 1000, dt = 0.0002;
const rng = makeRng(42);
const ts = [], meas = [];
let peak = 0;
for (let k = 0; k < N; k++) {
  const t = k * dt;
  const i = current(t, trueR, trueL, Um, w);
  ts.push(t); meas.push(i);
  peak = Math.max(peak, Math.abs(i));
}
for (let k = 0; k < N; k++) meas[k] += 0.01 * peak * gauss(rng);

let sseCalls = 0;
function sse(p) {
  sseCalls++;
  const [R, L] = p;
  if (R <= 0 || L <= 0) return 1e18;
  let s = 0;
  for (let k = 0; k < N; k++) {
    const d = current(ts[k], R, L, Um, w) - meas[k];
    s += d * d;
  }
  return s;
}

function nelderMead(fn, x0, step, iters = 300) {
  const n = x0.length;
  let simplex = [x0.slice()];
  for (let i = 0; i < n; i++) {
    const x = x0.slice(); x[i] += step[i]; simplex.push(x);
  }
  let vals = simplex.map(fn);
  for (let it = 0; it < iters; it++) {
    const order = vals.map((v, i) => i).sort((a, b) => vals[a] - vals[b]);
    simplex = order.map(i => simplex[i]); vals = order.map(i => vals[i]);
    const c = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) c[j] += simplex[i][j] / n;
    const worst = simplex[n];
    const pt = (a) => c.map((cj, j) => cj + a * (worst[j] - cj));
    const xr = pt(-1), fr = fn(xr);
    if (fr < vals[0]) {
      const xe = pt(-2), fe = fn(xe);
      if (fe < fr) { simplex[n] = xe; vals[n] = fe; } else { simplex[n] = xr; vals[n] = fr; }
    } else if (fr < vals[n - 1]) {
      simplex[n] = xr; vals[n] = fr;
    } else {
      const xc = pt(fr < vals[n] ? -0.5 : 0.5), fc = fn(xc);
      if (fc < Math.min(fr, vals[n])) { simplex[n] = xc; vals[n] = fc; }
      else {
        for (let i = 1; i <= n; i++) {
          simplex[i] = simplex[i].map((x, j) => simplex[0][j] + 0.5 * (x - simplex[0][j]));
          vals[i] = fn(simplex[i]);
        }
      }
    }
  }
  let b = 0;
  for (let i = 1; i <= n; i++) if (vals[i] < vals[b]) b = i;
  return { x: simplex[b], f: vals[b] };
}

const start = [3.0, 0.03];
const res = nelderMead(sse, start, [0.5, 0.005]);
const [Rfit, Lfit] = res.x;
const rmse = Math.sqrt(res.f / N);

console.log("\n=== РЕЗУЛЬТАТ ПІДБОРУ ===");
console.log(`Початкове наближення: R=${start[0]} Ом, L=${start[1]} Гн`);
console.log(`Істинні параметри:    R=${trueR} Ом, L=${trueL} Гн`);
console.log(`Знайдено:             R=${Rfit.toFixed(4)} Ом, L=${Lfit.toFixed(5)} Гн`);
console.log(`Похибка R: ${(100 * Math.abs(Rfit - trueR) / trueR).toFixed(2)} %, ` +
            `L: ${(100 * Math.abs(Lfit - trueL) / trueL).toFixed(2)} %`);
console.log(`СКВ залишків: ${rmse.toFixed(3)} А (шум ≈ ${(0.01 * peak).toFixed(3)} А)`);

console.log("\n=== ВИТРАТИ ===");
console.log(`Викликів цільової функції: ${sseCalls}, відліків N = ${N}, ` +
            `разом обчислень струму: ${sseCalls * N}`);
