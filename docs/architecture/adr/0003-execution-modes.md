# ADR 0003: Operation Execution Modes

- **Status**: Accepted
- **Date**: 2026-05-08

## Context
Different mesh operations require different visibility scopes (local chunk vs global statistics vs full topology).

## Decision
Standardize three execution modes:
1. **Streaming**: local chunk transforms with minimal memory.
2. **Two-pass**: pass-1 global analysis; pass-2 transform/write.
3. **Full-index**: build indexes/adjacency before execution.

Pipeline orchestration selects mode per operation contract.

## Alternatives Considered
1. Single full-in-memory mode: simplest implementation but poor memory scaling.
2. Streaming-only mode: memory efficient but insufficient for global/topology operations.
3. Dynamic mode inference without explicit operation declarations: flexible but less predictable and harder to validate.

## Consequences
### Positive
- Explicit capability modeling for extension authors.
- Predictable performance/memory profiles.
- Supports both simple and advanced operations.

### Negative
- More orchestration complexity in the pipeline.
- Requires disciplined operation metadata and tests.

## Follow-up
- Define operation manifest schema fields required for mode selection.
- Implement per-mode KPI instrumentation.
