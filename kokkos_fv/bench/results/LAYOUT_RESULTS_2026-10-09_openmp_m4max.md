# Layout benchmark results: Kokkos OpenMP backend, Apple M4 Max (2026-10-09)

**Raw data:** `layout_20261009-172317.csv` (240 rows).

**Tools:**
- `miniapp/fv_layout_bench`: standalone Kokkos 4.7.4 from the `moose-next` env, OpenMP,
  `OMP_PROC_BIND=spread OMP_PLACES=threads`.
- Meshes dumped from MOOSE `FaceInfo`/`ElemInfo` with the test object `FVMeshDumper`.

**Problem.** The same advection-diffusion-source system as `inputs/advdiff_kokkos.i`, with
coefficients evaluated inline (no parsed-function interpreter). Times are the median of 10
assemblies, covering the face/boundary passes plus the cell pass. Matrix zeroing and the PETSc
upload are not included.

**Correctness:**
- V0 (the replica of today's scheme) matches the matrix and RHS assembled by MOOSE's Kokkos path:
  exactly on orthogonal hex, and to about 1e-15 on tet and skewed hex.
- Every variant matches V0 to at most 6.6e-16 relative, in all 240 configurations.

## Variants

| ID | Layout / algorithm |
|---|---|
| V0 | **Today's scheme.** Padded per-(side, cell) geometry, one thread per directed side, coefficients evaluated per contribution, column found by binary search, atomics, separate matrix/RHS loops per kernel |
| V1 | Face SoA, one thread per face, coefficients once per face; still binary search + atomics |
| V2 | V1 + **precomputed CSR slots** (no search); atomics |
| V3 | Face pass → per-face scratch, off-diagonals into their unique slots; **per-cell gather** over `elem_faces` (no atomics) |
| V4 | Face pass → **half-face** `(side, cell)` scratch; gather writes the whole row (no atomics) |

## 1. Per-assembly time (ms) and speedup over V0, native ordering

| Mesh (cells) | Threads | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|---|
| hex 262k | 1 | 15.7 | 12.8 (1.23x) | 6.25 (2.51x) | **5.84 (2.69x)** | 9.05 (1.73x) |
| hex 262k | 8 | 2.43 | 2.48 (0.98x) | **1.67 (1.45x)** | 1.84 (1.32x) | 2.42 (1.00x) |
| hex 2.1M | 1 | 129.5 | 103.5 (1.25x) | 52.5 (2.47x) | **46.7 (2.77x)** | 51.6 (2.51x) |
| hex 2.1M | 4 | 33.7 | 27.1 (1.25x) | **13.6 (2.48x)** | 15.3 (2.21x) | 17.7 (1.90x) |
| hex 2.1M | 8 | 17.6 | 14.5 (1.21x) | **7.45 (2.36x)** | 11.9 (1.48x) | 13.8 (1.28x) |
| tet 1.5M | 1 | 60.8 | 62.5 (0.97x) | 34.7 (1.75x) | **33.4 (1.82x)** | 35.2 (1.72x) |
| tet 1.5M | 4 | 16.2 | 16.9 (0.96x) | **9.58 (1.69x)** | 11.9 (1.36x) | 10.9 (1.48x) |
| tet 1.5M | 8 | 9.51 | 9.54 (1.00x) | **8.62 (1.10x)** | 11.2 (0.85x) | 9.98 (0.95x) |

With RCM cell numbering and owner-sorted faces, the tet 1.5M case at 8 threads gives V2 4.36 ms
(2.2x over V0 native), V3 6.58 ms (1.4x) and V4 6.82 ms (1.4x).

## 2. Determinism: distinct matrix/RHS hashes over 10 assemblies

| | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| 1 thread, all meshes | 1 | 1 | 1 | 1 | 1 |
| 4-8 threads, skewed hex and tet | 5-10 | 10 | 10 | **1** | **1** |
| 4-8 threads, generated hex (native order) | 1 | 1 | 1 | 1 | 1 |

On structured, natively ordered hex meshes, the atomic variants happened to be reproducible:
threads rarely collide on the same entry. On the skewed hex and tet meshes they were not. The two
gather variants were bitwise reproducible everywhere.

## 3. Cell numbering and face order (8 threads, ms)

| Variant | hex 2.1M native/native | hex scrambled/native | hex rcm/owner | tet native/native | tet rcm/owner |
|---|---|---|---|---|---|
| V0 | 17.6 | 18.1 | 17.4 | 9.51 | 8.69 |
| V2 | 7.45 | **32.0** | 8.91 | 8.62 | **4.36** |
| V3 | 11.9 | **33.5** | 11.5 | 11.2 | **6.58** |
| V4 | 13.8 | **28.8** | 12.5 | 9.98 | **6.82** |

## 4. Memory for indexing/geometry structures (bytes per cell)

| Structure | hex | tet |
|---|---|---|
| Padded directed-side geometry + worklists (V0, today) | 528 | 352 |
| Face SoA (V1-V4) | 280 | 185 |
| `elem_faces` | 24 | 16 |
| Precomputed slots | 52 | 36 |
| Per-face scratch (V3) | 98 | 64 |
| Half-face scratch (V4) | 144 | 96 |

## 5. Measurement quality

- **1-thread results are stable.** The hex 2.1M and skewed-hex 2.1M meshes have identical
  connectivity, and they agree within about 2-5%.
- **Some 8-thread results are noisy.** Example: V0 skewed-hex morton/owner gave 40.3 ms against
  17.6 ms for the identical-connectivity hex case, and V0 should be insensitive to ordering.
- The M4 Max mixes performance and efficiency cores, and `spread` binding with 8 of 14 threads
  can place work on efficiency cores. Ordering conclusions at 8 threads are therefore drawn only
  where the effect is large (2x or more) and consistent across meshes.
- Each configuration was run once (10 timed assemblies). Repeats would tighten the 8-thread
  numbers.

## 6. Conclusions (CPU backend; GPU still to be measured)

1. **Precomputed matrix slots are the single biggest win.** V1 → V2 is 2x. Removing the binary
   search accounts for most of the gain over today's scheme: 2.4-2.5x on hex, 1.7-2.2x on tet.
   **Adopt (H3).**
2. **Compute-once face storage on its own gives little** (V1: 1.2x on hex, none on tet). With
   atomics, a face thread writes four matrix entries, so contention grows. It pays off together
   with slots, and it halves the geometry memory (528 → 280 B/cell on hex). **Adopt for memory
   and as the enabler of the rest (I1).**
3. **The gather (V3/V4) is the only deterministic option and remains faster than today's scheme.**
   On hex: 2.8x at 1 thread, 1.5x at 8 threads. But on the CPU it is slower than V2 at 8 threads,
   because the extra pass costs more than the CPU's cheap atomics. This matches the research
   (atomics are cheap on CPUs and recent GPUs). **Keep both modes (R3).** The deciding
   measurement has to be on a GPU.
4. **Half-face scratch (V4) did not beat per-face scratch (V3)** on the CPU, except on tet with
   RCM/owner where they are equal. **Hypothesis R2 is not supported on the CPU.** It needs to be
   re-tested on a GPU, where coalescing matters more than cache reuse.
5. **Ordering matters a lot for the new layouts, not for today's.** The fast variants are
   memory-bound. Scrambled numbering slows them 2.5-4.5x, while V0 barely changes because its
   cost is dominated by the search and recomputation.
   - libMesh's native order on generated meshes is already near-optimal. Morton or RCM gained
     nothing there.
   - **Renumbering must be combined with owner-sorted faces.** Morton or RCM with the native face
     order was often worse than no renumbering.
   - RCM + owner-sorted faces recovered native performance on hex and was best on tet (2x for
     V2). **Adopt RCM + owner sort (I6/I7) for unstructured and imported meshes.**
6. **Thread scaling from 1 to 8:**
   - V0 7.4x and V2 7.0x on hex;
   - V3 3.9x and V4 3.7x, because their extra pass is memory-bound.
   On a CPU, the gather's extra memory traffic is the limiting factor.

## 7. Not measured here

- **GPU behavior:** coalescing, atomic throughput, launch latency. This is the decisive test for
  atomics vs gather and for per-face vs half-face scratch. The same mini-app builds for CUDA
  unchanged (Kokkos views and lambdas).
- **The PETSc upload and owned-rows COO** (H4), and fence removal (H1). Those need the MOOSE
  integration.
- **The parsed-function cost** (the inline coefficients here isolate layout effects).
