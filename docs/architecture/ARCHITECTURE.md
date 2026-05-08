# Mesh Post Processor (`meshpp`) Architecture

## 0. Document Control

- **Project**: Mesh Post Processor (`meshpp`)
- **Version**: `v0.1-draft`
- **Status**: Draft
- **Authors**: TBD
- **Reviewers**: TBD
- **Last updated**: 2026-05-08
- **Related ADRs**:
  - `docs/architecture/adr/0001-language-toolchain.md`
  - `docs/architecture/adr/0002-data-layout.md`
  - `docs/architecture/adr/0003-execution-modes.md`

---

## 1. Purpose & Scope

### 1.1 Purpose
Define the architecture for a high-performance, extensible `.post.msh` post-processing tool that reads, transforms, and writes meshes with millions of nodes/elements.

### 1.2 In Scope
- `.post.msh` parsing, validation, and writing
- Operation pipeline and extension architecture
- Performance, memory, and reliability strategy
- CLI-driven execution model

### 1.3 Out of Scope (initially)
- GUI
- Distributed processing cluster
- Binary mesh format
- Remote storage backends

---

## 2. Stakeholders & Use Cases

### 2.1 Stakeholders
- Mesh processing developers
- Simulation engineers
- CI/CD and DevOps maintainers

### 2.2 Primary Use Cases
1. Apply coordinate transforms to large meshes.
2. Run chained operations from CLI or config.
3. Guarantee deterministic outputs for reproducible simulation workflows.
4. Add operations without parser/writer core changes.

---

## 3. Requirements

### 3.1 Functional Requirements
- FR-1: Read valid `.post.msh` files (header, comments, coordinates, elements).
- FR-2: Apply one or multiple operations in order.
- FR-3: Write valid `.post.msh` output.
- FR-4: Validate operation parameters and return actionable diagnostics.
- FR-5: Support execution modes (`streaming`, `two-pass`, `full-index`).

### 3.2 Non-Functional Requirements
- NFR-1: Handle meshes with millions of nodes/elements.
- NFR-2: Deterministic output for the same input + operation config.
- NFR-3: Memory usage bounded by configurable limits.
- NFR-4: Throughput-oriented buffered I/O.
- NFR-5: Testability and benchmarkability in CI.

### 3.3 Compatibility Constraints
- Text `.post.msh` compatibility with GiD post mesh conventions.
- Preserve IDs and ordering unless operation explicitly modifies them.

---

## 4. Format Model & Invariants

### 4.1 Logical File Structure
- Header: `MESH "..."`
- Optional metadata/comments (for example, `# color R G B`)
- `Coordinates ... End Coordinates`
- `Elements ... End Elements`

### 4.2 Core Invariants
- `dimension` matches coordinate arity.
- `Nnode` matches per-element connectivity arity.
- Every element node reference resolves to an existing node.
- Node IDs are unique within coordinates section.
- Element IDs are unique within elements section.

### 4.3 Error Taxonomy
- Syntax errors (malformed section/token)
- Semantic errors (invalid references/arity mismatch)
- Operation precondition errors
- Runtime resource errors (I/O failure, OOM)

---

## 5. Architecture Overview

### 5.1 Layered Modules
1. **`core/format`**: tokenizer, parser, serializer
2. **`core/model`**: mesh domain objects and data layout abstraction
3. **`core/io`**: buffered I/O and optional mmap adapters
4. **`core/ops`**: operation interface and registry
5. **`app/pipeline`**: stage orchestration and mode selection
6. **`cli`**: command-line interface

### 5.2 Dependency Rules
- `cli -> app/pipeline -> core/*`
- `core/*` must not depend on `cli`
- Operations depend on `core/model` contracts, not parser internals

---

## 6. Data Model

### 6.1 MeshHeader
- `name`
- `dimension`
- `elem_type`
- `nodes_per_element`

### 6.2 NodeTable
- Node IDs
- Coordinates (SoA preferred)
- Optional lazy ID->index map

### 6.3 ElementTable
- Element IDs
- Connectivity in contiguous flat storage
- Optional lazy adjacency/index structures

### 6.4 Metadata
- Preserved comments and recognized metadata entries

### 6.5 Data Layout Decision
- Default: SoA for coordinates, contiguous connectivity
- Rationale: cache locality and vectorization-friendly access patterns

---

## 7. Execution Model

### 7.1 Pipeline
`Read -> Validate -> Analyze -> Transform -> Write`

### 7.2 Operation Modes
- **Streaming**: chunk-local operations
- **Two-pass**: global analysis then transform
- **Full-index**: build indexes/adjacency before execution

### 7.3 Scheduling
- Parser: ordered
- Transform: parallel only for operations that declare thread safety
- Writer: ordered and deterministic

### 7.4 Determinism Rules
- Stable operation order
- Ordered chunk commit in writer
- Stable numeric formatting policy

---

## 8. Operation Extension Contract

### 8.1 Operation Manifest
- Name
- Version
- Parameter schema
- Required sections
- Required mode
- Thread-safety declaration

### 8.2 Lifecycle Hooks
- `configure(params)`
- `analyze(summary/context)` (optional)
- `execute(chunk/context)`
- `finalize()`

### 8.3 Preconditions/Postconditions
- Explicitly documented per operation
- Runtime-enforced by pipeline

### 8.4 Failure Behavior
- Recoverable warning vs fatal error policy
- Standardized error code mapping

---

## 9. I/O & Performance Strategy

### 9.1 Input
- Large buffered sequential reads
- Fast numeric parsing path

### 9.2 Output
- Large buffered writes
- Batched line emission
- Configurable numeric precision (v1 default policy: preserve in-memory `double` values and emit with stream default formatting; no fixed precision contract yet)

### 9.3 Large Mesh Handling
- Configurable chunk size
- Bounded stage queues
- Memory budget enforcement

### 9.4 Performance KPIs
- Throughput (MB/s)
- Time per million nodes/elements
- Peak RSS
- Scaling with thread count

---

## 10. Configuration & CLI

### 10.1 Initial Commands
- `meshpp apply --in <file> --out <file> --op <spec>`
- `meshpp apply --in <file> --out <file> --pipeline <yaml>`

### 10.2 Config Schema
- Ordered pipeline
- Global execution settings
- Operation-specific parameters

### 10.3 Validation Rules
- Strict schema validation before execution
- Actionable diagnostics on mismatch
- Metadata policy: parser is strict on required structural tokens (`MESH`, `Coordinates`, `Elements` blocks) and lenient on comments (`# ...`) by ignoring non-structural comment lines during parsing

---

## 11. Observability

### 11.1 Logging
- Levels: error/warn/info/debug/trace
- Structured fields: stage/op/chunk/time/memory

### 11.2 Metrics
- Parse/read/write durations
- Per-operation durations
- Peak memory
- Error counts by class

### 11.3 Progress Reporting
- Section-level and chunk-level progress
- Optional periodic status lines for long runs

---

## 12. Testing Strategy

### 12.1 Unit Tests
- Tokenizer/parser/serializer
- Data-model invariants
- Operation-level tests

### 12.2 Integration Tests
- End-to-end apply workflow
- Golden-file roundtrip
- Multi-operation chains

### 12.3 Fuzz/Property Tests
- Malformed-input parser resilience
- Connectivity and ID invariants

### 12.4 Benchmark/Regression
- Million-scale synthetic fixtures
- KPI threshold alerts in CI

---

## 13. Security & Robustness

- Defensive parsing and bounds checks
- Input size limits and caps
- Safe failure semantics for output files
- Sanitizers in CI (ASan/UBSan/TSan where applicable)

---

## 14. Deployment & CI/CD

- Toolchain matrix (Linux/macOS; GCC/Clang)
- Build profiles: Debug/RelWithDebInfo/Release
- CI stages:
  1. lint/static checks
  2. unit + integration
  3. sanitizer jobs
  4. benchmark smoke checks

---

## 15. Risks & Mitigations

- **Memory growth from adjacency-heavy operations** -> on-demand index construction and memory budget guards.
- **Parser bottlenecks** -> benchmark-driven optimization backlog.
- **Parallel nondeterminism** -> ordered writer and deterministic scheduling.

---

## 16. Roadmap Alignment

- M1: Foundation
- M2: Pipeline + basic operations
- M3: Scale and optimization
- M4: Extension ecosystem and contributor UX

---

## 17. Open Questions

- Compression support in v1 (`.gz`/`.zst`) or post-v1?
- (Closed) Numeric precision and formatting policy: preserve parsed `double` values and write with default stream formatting in v1.
- (Closed) Strict vs lenient policy: strict for structural sections/tokens, lenient for non-structural comments.
- Static operation registration only, or future dynamic loading?
