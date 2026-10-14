# Green-Gauss gradient variants: Kokkos OpenMP backend, Apple M4 Max (2026-10-14)

**Raw data:** `gradient_20261014-102129.csv` (432 rows).

**Setup:**
- Tool: `miniapp/fv_layout_bench --variants none --grad-variants G0,G1,G2,G3`, driven by
  `miniapp/run_gradient_bench.sh`.
- 20 timed passes per process; 3 processes per configuration; the tables show the median of the
  3 process medians.
- `OMP_PROC_BIND=spread OMP_PLACES=threads`.
- Meshes from MOOSE `FaceInfo`/`ElemInfo` dumps: hex 262k and 2.1M cells, skewed hex 2.1M (same
  connectivity as the 2.1M hex), tet 1.5M.

**Field:** `phi = 1 + xy + z^2 + sin(3x)` at cell centroids. Face values follow the CPU
`LinearFVGradientManager`:
- internal faces: linear interpolation with `FaceInfo::gC()`;
- Dirichlet faces: `phi` at the face centroid;
- other boundary faces: the cell value.

## Variants

All write the gradient with the components of a cell together (`3 * c + d`).

| ID | Algorithm | Inputs per cell (hex) |
|---|---|---|
| G0 | Per-cell gather over the padded per-(side, cell) geometry, as in today's Kokkos mesh | 432 B/cell of side data |
| G1 | Per-cell gather over the face SoA through `elem_faces`, the layout of the design goals | 256 B/cell face SoA + 24 B/cell `elem_faces` |
| G2 | One thread per face, atomic adds into both cells, then a volume pass | face SoA |
| G3 | G2 with 1/V folded into the atomic adds, a single face pass | face SoA |

## Correctness

- **Accuracy:** every variant matches the serial host reference to at most 8.4e-14 relative, on
  every mesh, ordering and thread count (exactly on the orthogonal hex).
- **Determinism:** G0 and G1 were bitwise reproducible everywhere.
- G2 and G3 with 4-8 threads on the skewed hex and the tet gave **13-20 distinct results out of
  20 passes**. On the generated hex they happened to repeat, as in the assembly study.

## Times (ms per gradient, median)

Skewed hex 2.1M has the same connectivity as hex 2.1M. It is used for the 8-thread comparison
because the hex 2.1M 8-thread runs were noisy (spread up to 52%).

| Mesh | Threads | Order | G0 | G1 | G2 | G3 |
|---|---|---|---|---|---|---|
| hex 262k | 1 | native | **1.81** | 2.45 | 4.34 | 4.72 |
| hex 262k | 8 | native | **0.52** | 0.60 | 1.03 | 1.11 |
| hex 2.1M | 1 | native | **14.6** | 20.2 | 35.3 | 38.3 |
| hex 2.1M | 4 | native | **4.63** | 5.82 | 9.35 | 9.97 |
| skewed hex 2.1M | 8 | native | **2.84** | 3.13 | 5.62 | 6.73 |
| skewed hex 2.1M | 8 | rcm/owner | **2.91** | 3.04 | 5.43 | 5.55 |
| tet 1.5M | 1 | native | **7.66** | 10.7 | 19.4 | 20.4 |
| tet 1.5M | 8 | native | **1.57** | 2.11 | 3.12 | 3.05 |
| tet 1.5M | 8 | rcm/owner | **1.63** | 1.73 | 2.92 | 2.91 |
| hex 2.1M | 1 | scrambled | **32.9** | 233 | 127 | 142 |
| tet 1.5M | 8 | scrambled | **3.15** | 13.9 | 10.4 | 10.8 |

## Findings (CPU only)

1. **Gathers beat atomics by 1.5-1.8x** at every thread count and native or RCM ordering, and
   they are the only deterministic option.
   - The per-face atomic variants do 6 atomic read-modify-writes per face, against 3 plain stores
     per cell for a gather.
   - Folding 1/V into the atomics (G3) is slower than a separate volume pass (G2): two divides
     per face cost more than one more pass over the cells.
2. **G0 (padded side data) is the fastest on the CPU: 1.0-1.4x faster than G1 (face SoA).**
   - G0 streams a cell's side data, contiguous because the padded arrays are side-fastest. It
     reads only `phi(c)` and `phi(neighbor)`.
   - G1 goes through `elem_faces` to the face arrays, a second indirection, and reads `phi` at
     both face cells.
   - The cost of G0 is memory: 432 against 280 B/cell on hex.
3. **G1 is the most sensitive to cell numbering.**
   - With scrambled numbering G1 is 7-12x slower than with native numbering. The face data of a
     cell is then spread across memory, and G1 even falls behind the atomic variants.
   - G0 only loses about 2x, because its geometry stays contiguous per cell.
   - RCM with owner-sorted faces closes most of the G1-G0 gap on the tet at 8 threads (1.73
     against 1.63 ms).
   - Generated meshes are already well numbered; imported meshes need RCM plus owner sorting
     (I6/I7) for G1.
4. **Cost relative to assembly:** about 3 ms per gradient at 2.1M cells and 8 threads. That is
   about 40% of the mini-app's best assembly (V2, 7.45 ms) and about 2% of MOOSE's current Kokkos
   assembly (134 ms). The gradient pass is not a bottleneck on the CPU.

## What the GPU must decide

- **Coalescing.** On a GPU, consecutive threads (cells) read the padded side arrays with a stride
  of the side count (6 for hex). The side-fastest layout that helps the CPU stream a cell's data
  may coalesce poorly; a cell-fastest layout of the same data would coalesce. G1's face arrays
  are read through an index, so coalescing depends on the numbering.
- **Atomics.** Whether atomics are still 1.5-1.8x behind on the GPU, given native fp64 atomics in
  L2.
- **Memory.** G0's extra 150 B/cell matters more on a GPU, where memory limits the cells per
  device.

**Suggested GPU additions:**
- a G0 variant with cell-fastest padded arrays;
- a G1 variant with the face arrays of each cell copied into a per-cell (CSR-ordered) structure,
  which trades G0's memory for G1's single storage.
