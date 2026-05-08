# Adding a new meshpp operation

This guide explains how to add a new operation to the `meshpp` pipeline.

## 1) Implement the operation class

Add a class in `MESHPP/src/operations.cc` deriving from `MeshOperation`:
- implement `Name()` with a unique operation name
- implement `Configure(spec)` to parse and validate arguments
- implement `Apply(MeshData*)` to transform mesh data

Use `ValidationReport` with `ExitCode::kUsageError` for invalid operation specs.

## 2) Register the operation

Update `CreateOperation(...)` in `MESHPP/src/operations.cc` to construct your class when the operation name matches.

## 3) Keep behavior deterministic

Operations are applied in user-provided order. Avoid non-deterministic iteration or RNG unless explicitly seeded and documented.

## 4) Add tests

- Positive pipeline test in `MESHPP/tests/run_apply_pipeline.sh` style.
- Negative-case test in `MESHPP/tests/run_negative_cases.sh` for invalid specs or malformed input.

## 5) Validate performance impact

Run:

```bash
cd MESHPP
./benchmarks/run_benchmark.sh 10000
```

Check stage timings (`stats.*`) for regressions.
