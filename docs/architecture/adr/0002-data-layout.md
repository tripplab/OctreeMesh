# ADR 0002: Core Mesh Data Layout

- **Status**: Accepted
- **Date**: 2026-05-08

## Context
The system must process millions of nodes/elements efficiently and support both streaming and indexed operations.

## Decision
Use a data model centered on:
- **SoA coordinates** for nodes (`x[]`, `y[]`, `z[]` by dimension)
- **Flat contiguous connectivity** for elements
- Optional **lazy** ID->index maps and adjacency/index structures

## Alternatives Considered
1. AoS node structs (`struct Node {id,x,y,z;}`): easier ergonomics but weaker cache behavior for coordinate-heavy loops.
2. Nested vectors for connectivity: simpler but fragmented memory and weaker locality.
3. Always-on global indexes: convenient but excessive memory overhead for streaming operations.

## Consequences
### Positive
- Better cache locality and SIMD-friendly traversal.
- Reduced memory fragmentation.
- Flexibility for both streaming and topology-heavy operations.

### Negative
- Slightly higher complexity in abstractions and API design.
- Requires careful validation around index/map coherence.

## Follow-up
- Define precise ownership/lifetime rules for optional indexes.
- Add benchmarks comparing streaming vs indexed operation paths.
