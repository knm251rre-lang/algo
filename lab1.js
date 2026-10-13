function rectangle(a, b) {
  if (!(a > 0) || !(b > 0)) throw new Error("Сторони мають бути додатними");
  const area = a * b;
  const perimeter = 2 * (a + b);
  return { area, perimeter };
}

function gcd(a, b) {
  a = Math.abs(a);
  b = Math.abs(b);
  while (b !== 0) {
    const r = a % b;
    a = b;
    b = r;
  }
  return a;
}

function gcdTrace(a, b) {
  a = Math.abs(a);
  b = Math.abs(b);
  let step = 1;
  while (b !== 0) {
    const r = a % b;
    console.log(`  крок ${step++}: ${a} mod ${b} = ${r}`);
    a = b;
    b = r;
  }
  return a;
}

console.log("=== Прямокутник ===");
for (const [a, b] of [[3, 4], [5.5, 2], [10, 10]]) {
  const r = rectangle(a, b);
  console.log(`a=${a}, b=${b} -> S=${r.area}, P=${r.perimeter}`);
}

console.log("\n=== Алгоритм Евкліда ===");
for (const [a, b] of [[48, 18], [1071, 462], [17, 5], [100, 0]]) {
  console.log(`НСД(${a}, ${b}):`);
  const g = gcdTrace(a, b);
  console.log(`  результат = ${g}`);
}
