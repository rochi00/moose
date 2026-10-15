# Phase 0 implementation plan: face-based Kokkos linear FV assembly

**Started:** 2026-10-14. **Design:** `kokkos_ns_linearfv_detailed_design.md` sections 3-5 and 14.

Phase 0 moves the existing Kokkos linear FV objects from one thread per directed element side,
reading the side-fastest per-(side, element) arrays of the shared Kokkos mesh, to one thread per
face over a face array owned by the linear FV layer. The GPU campaign measured the target: the
atomic face pass (V2, fused V2F) was 1.4-2.1x faster than today's scheme on an A100.

## Decisions (2026-10-14)

| Question | Decision |
|---|---|
| Where the face and cell arrays live | A new `Moose::Kokkos::FVMesh`, owned by the linear FV layer and built from `FaceInfo`/`ElemInfo`. The shared Kokkos `Mesh` and FE objects are unchanged. |
| Today's per-directed-side hooks | Replaced. All Kokkos linear FV objects move to the per-face hooks of design 5.2; the old hooks are removed. |
| Assembly mode | Atomic (default) and gather (opt-in, bitwise reproducible), both in phase 0. |

## Steps

| Step | Content | Done when |
|---|---|---|
| 1. `FVMesh` | Face arrays: two cell indices, the side index on each, area times coordinate factor, unit normal, `gC`, `\|dCN\|`, centroid, boundary ID. Cell arrays: volume times coordinate factor, centroid, cell-fastest `elem_faces`. Ghost cells use the Kokkos mesh's contiguous IDs. Built in `FEProblemBase::initKokkos()` and after mesh changes. | A test object compares every value with `FaceInfo`/`ElemInfo`, serial and with `-p 3` |
| 2. Face lists and slots | Per flux kernel: its internal faces (variable active on both sides) and boundary faces (with a Kokkos BC), and the COO slots of the 2x2 face block. Rows are written only for owned cells. | Slots match `getMatrixSlot()` (debug check) |
| 3. Hooks and datums | `FaceDatum`, `CellDatum`; `computeInternalFaceContribution` returning `{elem, neighbor, rhs_elem, rhs_neighbor}`, `computeBoundaryFaceContribution` returning `{diag, rhs}`, `computeCellContribution` returning `{diag, rhs}`. Face coefficients stored per face. Migrate diffusion, advection, source and the implicit Euler test object. | Kokkos linear FV tests pass with their current gold files |
| 4. BCs | BC relations computed over the BC's faces from `FVMesh` | BC tests pass |
| 5. Gather mode | Per-face scratch per flux kernel and a cell gather in local side order; `Problem/kokkos_linear_fv_assembly = atomic|gather` | Gather results reproducible across runs and thread counts; tests pass in both modes |
| 6. Measure | CPU: parity, sweep, event counts. GPU: `overhead` and `miniapp` sections | Assembly time compared with the mini-app's V2 |

**Added on 2026-10-15** (first deferred): face ranges, BCs on `FVMesh` faces, and an
`FVSystemAssembly` that owns the matrix and right-hand side path. **Still deferred:** the G0C
gradient side table (phase 1), renumbering (not needed per the GPU results).

## Progress

- **Steps 1-5 done (2026-10-14).**
  - `Moose::Kokkos::FVMesh` (`KokkosFVMesh.h/.K`), owned by `FEProblemBase`, created when the
    first Kokkos linear FV object is added, built in `initKokkos()` and rebuilt after mesh changes.
    Test: `test/tests/kokkos/fvmesh` (`KokkosFVMeshTest`), serial, `-p 3` and distributed mesh.
  - `FaceDatum`, `CellDatum` (`KokkosFVDatum.h`). One entry point per object,
    `computeContributions(matrix, rhs)`; flux kernels loop over internal faces (four precomputed
    slots each) and boundary faces; elemental kernels over cells. The face-coefficient pass is gone:
    each face is evaluated once per assembly anyway. `FVDatum` is removed.
  - BCs keep their worklists and map each `(element, side)` to its face through `elem_faces`.
  - Gather mode: `Problem/kokkos_linear_fv_assembly = gather`.
- **Results (CPU, M4 Max):**
  - Kokkos linear FV tests pass with their gold files, serial, `-p 3`, distributed mesh, and in
    gather mode. `bench.py parity` is unchanged case by case.
  - Launches per assembly 16 to 12; fences unchanged at 54 (all in shared code).
  - Steady assembly 1.4-1.7x faster at 32k-2.1M cells and 1-8 threads; 2.1M cells on 8 threads
    134 to 82 ms (`sweep_20261014-162518.csv`).
  - Gather mode: 1 distinct matrix in 5 runs on 8 threads (tet, skewed), against 5 in atomic mode.
    It costs 1.2-1.3x the atomic assembly on the CPU (2.1M cells, 8 threads: 106 against 82 ms,
    `sweep_20261014-163936.csv`).
- **Not changed, by decision (2026-10-15):** the Kokkos mesh still builds its per-side geometry
  for linear FV problems (about 760 B per hex cell), although no linear FV object reads it now.
  `setNeedsElementSideGeometry()` also enumerates the one-layer ghost elements that `FVMesh` and
  the Kokkos system need, so dropping the geometry means changing the shared Kokkos mesh, which
  the Kokkos FE objects use. That change is left out.
- **Ranges, BCs and `FVSystemAssembly` done (2026-10-15).**
  - `FVMesh` orders its faces into ranges (`FaceRange`) keyed by elem block, kind (interior,
    processor, interface, boundary), neighbor block, which cells are owned, and the boundary IDs
    of each side; within a range by owned cell. Kernels and BCs keep lists of ranges
    (`FaceRangeList`); a thread finds its range by binary search over the range offsets. Per face,
    kernels store only matrix slots; the flip, the BC and its offset are per range.
  - BCs build their faces from the ranges where their variable lives on one side, the cell on
    that side is owned, and that side carries exactly one boundary ID, one of the BC's. Boundary
    IDs are kept per side: a sideset holds element sides, and `FaceInfo::boundaryIDs()` merges both
    sides (the block-restricted internal Neumann test has a different sideset on each side of the
    interface). This replaces the BC worklists of the shared boundary restriction.
  - `FVSystemAssembly` (`KokkosFVSystemAssembly.h/.K`), owned by `LinearSystem`: one value array
    per matrix tag laid out as the owned rows of the Kokkos system's sparsity, so the precomputed
    slots stay valid; the pattern is preallocated with the owned rows only; right-hand sides are
    added straight into the PETSc vector arrays (staged on device only for host vectors on a GPU).
    Objects get the device arrays of their tags before each dispatch, so the per-assembly tag flag
    copies are gone. The Kokkos system now only provides the solution vectors.
  - Events per assembly (CPU): fences 54 to 20, allocations 3 to 1, deep copies 4 to 0, launches
    12 to 10. The remaining fences are on the solution-vector path of the shared Kokkos system,
    except `KokkosFVWait` and the one after the upload.
  - Tests pass serial, `-p 3`, distributed mesh, and in gather mode; parity unchanged.
  - Steady assembly 1.05-1.17x faster than the face-based version alone, and 1.6-1.9x faster than
    before phase 0; 2.1M cells on 8 threads: 134 to 82 to 74 ms (`sweep_20261015-084411.csv`).
- Step 6 (GPU measurement) is next.
