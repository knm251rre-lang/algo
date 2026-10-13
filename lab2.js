function fibonacci(n) {
  if (n < 1) return -1;
  if (n === 1 || n === 2) return 1;
  let prev = 1, cur = 1;
  for (let i = 3; i <= n; i++) {
    const next = prev + cur;
    prev = cur;
    cur = next;
  }
  return cur;
}

const DIGITS = "0123456789ABCDEF";
function toBase(x, base) {
  if (!Number.isInteger(x) || x < 0) throw new Error("Потрібне невід'ємне ціле");
  if (x === 0) return "0";
  let result = "";
  while (x > 0) {
    result = DIGITS[x % base] + result;
    x = Math.floor(x / base);
  }
  return result;
}
const toBinary = (x) => toBase(x, 2);
const toHex = (x) => toBase(x, 16);

console.log("=== Фібоначчі (ітеративно) ===");
let line = [];
for (let n = 1; n <= 15; n++) line.push(fibonacci(n));
console.log(line.join(" "));
console.log("F(50) =", fibonacci(50));

console.log("\n=== Системи числення ===");
for (const x of [0, 5, 10, 255, 1024, 48879]) {
  console.log(`${x} -> bin ${toBinary(x)} | hex ${toHex(x)} | перевірка: ${x.toString(2) === toBinary(x) && x.toString(16).toUpperCase() === toHex(x)}`);
}
