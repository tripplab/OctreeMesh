#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
GEN="$ROOT_DIR/MESHPP/benchmarks/generate_hex_grid.py"
BIN="$ROOT_DIR/MESHPP/meshpp_apply"
TMP_DIR="${TMPDIR:-/tmp}/meshpp_bench"
mkdir -p "$TMP_DIR"

N_HEX="${1:-10000}"
IN="$TMP_DIR/in_${N_HEX}.post.msh"
OUT="$TMP_DIR/out_${N_HEX}.post.msh"

python3 "$GEN" "$N_HEX" "$IN"
if command -v /usr/bin/time >/dev/null 2>&1; then
  /usr/bin/time -f "wall_s=%e rss_kb=%M" "$BIN" --in "$IN" --out "$OUT" --op scale:1.01 --op translate:0.1,0.2,0.3 --perf_stats
else
  time "$BIN" --in "$IN" --out "$OUT" --op scale:1.01 --op translate:0.1,0.2,0.3 --perf_stats
fi
