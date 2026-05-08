#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

cd "$ROOT_DIR/MESH2STL"
./tests/ci.sh

cd "$ROOT_DIR/MESHPP"
./tests/ci.sh
