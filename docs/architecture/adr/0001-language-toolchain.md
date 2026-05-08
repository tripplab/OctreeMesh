# ADR 0001: Language & Toolchain Baseline

- **Status**: Accepted
- **Date**: 2026-05-08

## Context
The tool must process very large text mesh files with deterministic behavior and high throughput while remaining maintainable and extensible.

## Decision
Adopt **C++20** as baseline language standard, with **CMake** as build system and a CI matrix using GCC and Clang.

## Alternatives Considered
1. C++17: stable, but misses modern standard features that improve expressiveness and safety.
2. Rust: strong safety guarantees, but higher migration/onboarding overhead for the current code ecosystem.
3. Python + native extensions: faster iteration but weaker baseline throughput and more packaging complexity.

## Consequences
### Positive
- Modern language facilities for robust architecture.
- Strong performance potential for million-scale workloads.
- Good portability across Linux/macOS.

### Negative
- Increased implementation complexity vs scripting-language prototypes.
- Requires disciplined performance profiling and ABI/toolchain hygiene.

## Follow-up
- Define CI compiler versions.
- Establish lint/format/static-analysis policy.
