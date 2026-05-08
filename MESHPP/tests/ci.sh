#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

make clean
make
./tests/run_roundtrip.sh
./tests/run_apply_pipeline.sh
./tests/run_negative_cases.sh
