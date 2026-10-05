// Wider display bins average the exact 1-bp input; source values stay untouched.
export function displayedCoverage(row, scale, binSize, normalized) {
  if (!row || row.bin_size !== 1) return [];
  const source = normalized ? row.normalized_values : row.values;
  const result = [];
  for (let pos = scale.start; pos <= scale.end; pos += binSize) {
    const hiPos = Math.min(scale.end, pos + binSize - 1);
    let sum = 0;
    for (let i = pos - row.start; i <= hiPos - row.start; i++) sum += source[i] || 0;
    result.push({ position: pos, end: hiPos, value: sum / (hiPos - pos + 1) });
  }
  return result;
}
