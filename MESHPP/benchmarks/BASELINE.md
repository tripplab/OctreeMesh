# meshpp benchmark baseline (PR4/PR5)

Environment: local CI sandbox, `g++ -O2`, single process.

## Command

```bash
cd MESHPP
./benchmarks/run_benchmark.sh <N_HEX>
```

## Baseline captures

### N_HEX=10000

- `stats.read_ms=129.359`
- `stats.validate_ms=0.643887`
- `stats.operations_ms=0.901384`
- `stats.write_ms=89.1096`
- `stats.nodes=80000`
- `stats.elements=10000`
- `real=0.234s`

### N_HEX=50000

- `stats.read_ms=700.686`
- `stats.validate_ms=2.11121`
- `stats.operations_ms=7.42596`
- `stats.write_ms=474.863`
- `stats.nodes=400000`
- `stats.elements=50000`
- `real=1.215s`

## Usage note

Use these as regression anchors, not absolute thresholds across machines.
