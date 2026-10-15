# Kokkos Navier-Stokes Linear FV: Detailed Design

**Status:** proposal, partly implemented. On the branch `rc/kokkos_moose_navier_stokes`, the
existing Kokkos linear FV objects have:

- one thread per face over an `FVMesh` face array built from `FaceInfo` (section 3, without the
  face ranges), with the per-face hooks of section 5 and atomic or gather assembly (4.4);
- precomputed matrix slots;
- the empty host loop skipped, and an error on systems that mix Kokkos and host linear FV objects
  (section 10);
- one fence per assembly, with no device allocations per dispatch;
- device PETSc types by default on GPU builds;
- `--kokkos-threads`;
- only the assembled system synchronized per assembly.

See `findings_2026-10-12.md`, `findings_2026-10-14.md` and `phase0_implementation_plan.md`.
Nothing else here is implemented.

**Baseline:** MOOSE `next` at `410cc90c7e`. Paths are relative to the `moose/` repository root.

**Inputs:**
- `kokkos_ns_linearfv_design.md`: the overview design (decisions D1-D11, open questions Q1-Q6).
- `gpu_fv_design_research.md`: evidence from other GPU codes and revisions R1-R12.
- `ns_linearfv_design.md`: the CPU design being ported.
- `kokkos_moose_design_choices.md`: the Kokkos-MOOSE conventions this design must follow.

This document refines the overview into concrete data structures, algorithms, interfaces,
synchronization rules, and a verification plan. Section numbers here are independent of the
overview document.

---

## 1. Scope, targets, and non-goals

### 1.1 Scope

- A GPU-resident segregated SIMPLE/PIMPLE solver for weakly compressible and incompressible
  Navier-Stokes using linear FV and Rhie-Chow interpolation, with the same discretization as the
  CPU path in `modules/navier_stokes`.
- Framework support (mesh, system, kernel and BC infrastructure) that also serves the existing
  Kokkos linear FV objects.

### 1.2 Performance targets

| Target | Value | Source |
|---|---|---|
| Problem size per GPU where we expect good efficiency | 2-10M cells | SPUMA, NekRS, OGL, commercial codes |
| Problem size per GPU where we require correctness but expect poor efficiency | 0.1-1M cells | same |
| MPI ranks per GPU | 1 | OGL oversubscription result |
| Host-device synchronizations per SIMPLE iteration, excluding inside KSP | O(number of equations) scalar reads only | SPUMA sync overhead |
| Host-device copies of mesh-sized data inside the iteration loop | 0 | EXASIM, AMD offload |
| Halo exchanges per SIMPLE iteration | at most about 2 per equation plus 3 for the pressure-velocity coupling | PetscSF latency budget |
| Atomics in assembly | 0 (default mode) | determinism, miniAero |

### 1.3 Non-goals for this design

- Mesh adaptivity and displaced meshes (Kokkos-MOOSE rejects both today).
- Mixed precision.
- A block-coupled pressure-velocity solver. The data structures must not rule it out (section
  13.3), but it is not designed here.
- Porting fluid-property user objects. Properties are constants, parsed functions, or Kokkos
  cell-field expressions in the first phases.

---

## 2. Architecture overview

```
 Input file ([LinearFVKernels], [LinearFVBCs], [UserObjects], [FVProperties]*, Executioner)
     |
     v
 Kokkos objects (Moose::Kokkos namespace, registered with registerKokkos* macros)
     |  flux kernels, elemental kernels, BCs, property objects, KokkosRhieChowMassFlux
     v
 Kokkos FV system layer (per LinearSystem)          Kokkos FV mesh layer (per mesh)
     |  owned-row sparsity, COO slot maps,              |  cell SoA, face SoA, elem_faces,
     |  half-face scratch, cell scratch,                |  face ranges, renumbering,
     |  gather + fused relaxation, upload               |  half-face slot map
     v                                                  v
 PETSc device objects: Mat (aijkokkos/hypre), Vec (kokkos), KSP, PetscSF ghost updates
```

`[FVProperties]*` stands for whatever syntax section 7 settles on. It may simply be
`[FunctorMaterials]` with Kokkos objects.

**Data residency rule.** Everything indexed by cell, face, or nonzero lives on the device for the
whole run. The host holds only:

- range tables (a few hundred entries),
- scalars (time, dt, relaxation factors, norms),
- the host mirror of the mesh used during setup.

---

## 3. Mesh layer

### 3.1 Index spaces

| Index | Type | Range | Notes |
|---|---|---|---|
| `CellID c` | `ContiguousElementID` | `[0, n_local)` owned; `[n_local, n_local + n_ghost)` ghost | Owned cells are grouped by block and renumbered for locality (3.2). Ghost cells are the one-layer FV ghost set that already exists. |
| `FaceID f` | `MOOSE_KOKKOS_INDEX_TYPE` | `[0, n_faces)` | Every face with at least one owned cell. |
| `HalfFaceID h` | `MOOSE_KOKKOS_INDEX_TYPE` | `side + max_sides * c` for owned `c` | An (owned cell, side) pair. Indexes the `(side, cell)` LEFT-layout scratch. |
| `SlotID k` | `PetscCount` | `[0, nnz_owned)` | Position in the COO value array. |
| `RangeID r` | `unsigned int` | small | Index into the face range table. |

The existing `ContiguousElementID` stays the cell index, so FE and FV Kokkos code keep sharing
`_elem_info`, dof maps and block ranges.

### 3.2 Cell renumbering (R4)

The current numbering follows libMesh iteration order within each block (`KokkosMesh.K:36-140`).
The proposal adds an optional locality pass:

1. For each block, collect its owned cells and their intra-block, intra-rank face adjacency from
   libMesh neighbor pointers.
2. Order them with reverse Cuthill-McKee, or with a Morton key of the centroid. The Morton key is
   cheaper and parallel-friendly; the LBM study found Hilbert no better.
3. Assign contiguous IDs in that order, still block by block, so `subdomain_elem_id_ranges` keeps
   its meaning.

Notes:
- Ghost cells keep the current scheme (appended after owned cells).
- The renumbering is stored in the existing `kokkos_contiguous_elem_id` extra integer, so every
  consumer (FE Kokkos objects, dof maps, materials) sees it automatically.
- It changes the order of the Kokkos FE loops too. That is harmless for correctness but can
  shift floating-point summation order in FE atomics, so it must be guarded by a mesh parameter
  (for example `kokkos_cell_ordering = native|rcm|morton`, default `native` until benchmarked).

### 3.3 Cell SoA

```cpp
struct FVCellData                         // all indexed by CellID unless noted
{
  Array<ElementInfo> info;                // existing: type, id, contiguous subdomain
  Array<Real>        volume;              // V * coord factor; owned + ghost
  Array<Real3>       centroid;            // owned + ghost (ghost needed for RZ and reconstructions)
  Array2D<FaceID>    faces;               // (side, c), LEFT layout; invalid past n_sides; owned only
  Array2D<PetscCount> slot;               // (side, c): COO slot of the off-diagonal entry for that side,
                                          // or -1 if the side has no in-system neighbor; owned only
  Array<PetscCount>  diag_slot;           // COO slot of the diagonal; owned only
};
```

Ghost cell entries exist only where a face needs them: volume, centroid, and the dof index from
the existing system maps. Ghost cells get no `faces` or `slot` entries.

### 3.4 Face SoA

```cpp
struct FVFaceData                         // all indexed by FaceID
{
  Array<CellID>  elem;                    // canonical "elem" side (always valid)
  Array<CellID>  neighbor;                // invalid on domain boundaries; may be a ghost CellID
  Array<uint8_t> elem_side;               // libMesh side index on elem
  Array<uint8_t> neighbor_side;           // libMesh side index on neighbor (if any)
  Array<Real>    area;                    // A * face coord factor
  Array<Real3>   normal;                  // unit normal, outward from elem
  Array<Real>    gc;                      // FaceInfo::gC(); 1 on boundaries
  Array<Real3>   dcn;                     // C_N - C_E, or x_f - C_E on boundaries
  Array<Real>    inv_dn;                  // 1/|dCN . n| (NOC on) or 1/|dCN| (off); per-mesh flag
  // optional, allocated on request:
  Array<Real3>   noc;                     // n - eCN / (n . eCN)
  Array<Real3>   centroid;                // face centroid
  Array<Real3>   skew;                    // skewness-correction vector
};
```

**Where face data comes from.** All geometric values are copied from MOOSE `FaceInfo` objects so
that GPU results are comparable to the CPU (Q2 in the overview).

**The NOC flag.** `inv_dn` depends on whether non-orthogonal correction is used. The CPU allows
this per kernel. Allocate one array per distinct setting that is requested, so in practice at most
two arrays.

**Half-face mapping.** A face writes into its two half-faces:

- the elem side: `h_E = elem_side[f] + max_sides * elem[f]`;
- the neighbor side: `h_N = neighbor_side[f] + max_sides * neighbor[f]`, used only if the neighbor
  is owned.

Both indices are computed from stored data; no extra array is needed.

### 3.5 Face set construction

The faces are built on the host, once per mesh change, by a `Moose::Kokkos::FVMesh` that the
linear FV layer owns (decided 2026-10-14). The shared Kokkos `Mesh` and the FE objects that read it
are unchanged. `FEProblemBase` builds the `FVMesh` in `initKokkos()` when Kokkos linear FV objects
exist, after `setupFiniteVolumeMeshData()`, and rebuilds it after a mesh change once the `FaceInfo`
objects are rebuilt.

1. **Collect.** For each `FaceInfo` in `MooseMesh::faceInfo()`, keep the face if its elem or its
   neighbor is owned. This is the same set as the CPU `_face_info`.
2. **Orient.** Keep `FaceInfo`'s elem/neighbor assignment. It follows `Moose::FV::elemHasFaceInfo`,
   which uses global IDs, so it is canonical across ranks (D5). This requires no work.
3. **Map.** Convert elem and neighbor to `CellID` through the contiguous-ID extra integer for owned
   cells and `ghost_elem_id_mapping` for ghost cells. Every neighbor of an owned cell must resolve;
   otherwise error out (the ghosting is insufficient).
4. **Classify** into a `FaceKind`, with keys `(block_E, block_N)` or `(block_E, boundary)`:
   - `INTERIOR`: both cells owned, same block;
   - `PROCESSOR`: one cell is a ghost, same block;
   - `INTERFACE`: different blocks, both cells exist;
   - `BOUNDARY`: no neighbor.

   A boundary face with more than one boundary ID is an error, as on the CPU. Interior faces that
   carry a sideset are classified by cells as usual; a separate map from sideset to face list
   covers internal sidesets.
5. **Sort** by `(primary block, kind, secondary key, owned cell ID of the first owned side)`.

   The primary block of an interface or processor face is the block of its first owned side, so
   each rank's sort key depends only on data it has.
6. **Build range tables.** Emit `FaceRange{begin, end, kind, block, other}` for each contiguous
   group, plus a map from boundary ID to its list of ranges (one per adjacent block).
7. **Fill `faces(side, c)`** for owned cells, and check the invariants:
   - every owned cell side appears exactly once;
   - every interior face appears in two half-faces when both sides are owned.
8. **Fill the geometry arrays**, then copy everything to the device.

**Internal sidesets** (sidesets on interior faces, used by CHT and baffles) are a list per boundary
ID of interior `FaceID`s. BC objects on internal sidesets loop over this list rather than over a
range. Faces where the variable lives on only one side are handled through variable block
restriction (3.6), not through this list.

**As implemented (2026-10-15).** `FVMesh::update()` sorts by `(elem block, kind, neighbor block,
elem owned, neighbor owned, elem-side boundary IDs, neighbor-side boundary IDs)`, then by owned
cell. Splitting processor faces by which cell is owned, and interface faces by orientation, makes
the rows a kernel writes and the side a variable lives on uniform within a range. Boundary IDs are
kept per side, not as `FaceInfo::boundaryIDs()`: a sideset holds element sides, and the two sides
of an interface can carry different sidesets. A face with several boundary IDs on the side where a
variable lives is an error only for that variable, as on the CPU. Internal sidesets on faces with
the variable on both sides are ignored, also as on the CPU, so no separate list is built.

### 3.6 Variable views over the ranges

A variable is block-restricted to a set of blocks `B`. For a given variable, each range falls into
one category:

| Range kind | Variable on E's block | Variable on N's block | Treated as |
|---|---|---|---|
| INTERIOR / PROCESSOR | yes | (same block) | internal face |
| INTERFACE | yes | yes | internal face |
| INTERFACE | yes | no | boundary face on the E side |
| INTERFACE | no | yes | boundary face on the N side (flip the normal sign) |
| BOUNDARY | yes | n/a | boundary face |

At `initialSetup`, each flux kernel and BC builds a **range list** from this table (host and
device). This replaces the per-kernel `_internal_face_ids` and `_boundary_face_ids` vectors that
exist today.

Each range also carries `orientation = +1/-1`, so kernels evaluated on N-side boundaries use the
outward normal of the side where the variable lives. This replaces the CPU `ELEM`/`NEIGHBOR`
`faceType` branching.

**Thread mapping.** A launch over a range list of length `R` uses `n = sum(sizes)` threads. Thread
`t` finds its range by a branch-free binary search over a device prefix array of length `R + 1`
(`R` is typically 1-10), then computes `f = range.begin + (t - prefix[r])`.

### 3.7 Memory budget

The table below is for a hexahedral mesh: about 3 faces per cell, 6 sides, 7 nonzeros per row,
8-byte reals, 4-byte indices.

| Data | Bytes per cell | Notes |
|---|---|---|
| Face SoA, required fields | about 3 x 70 = 210 | elem, neighbor, sides, area, normal, gc, dcn, inv_dn |
| Face SoA, optional fields | up to 3 x 72 = 216 | noc, centroid, skew |
| Cell SoA | about 6 x 4 + 6 x 8 + 8 + 8 + 24 = 112 | faces, slot, diag_slot, volume, centroid |
| Gradient side table (6.2) | about 6 x 72 = 432 | neighbor, boundary ID, area, gc, normal, centroid per (side, cell); geometry only, shared by all variables |
| Half-face scratch per assembled system | 6 x 3 x 8 = 144 | diag, off, rhs per half-face |
| Cell scratch per assembled system | 16 | elemental diag and rhs |
| COO values per system (our array) | 7 x 8 = 56 | removed if we write into PETSc's array (13.2) |
| PETSc matrix per system (AIJ values + columns + COO maps) | 7 x (8 + 4 + about 16) = about 200 | PETSc-owned |
| Solution + 2 old states per variable | 24 | |
| Gradient, one state, per variable | 24 | |

**Scratch is shared.** Systems are assembled one after another, so the half-face and cell scratch
can be **one shared buffer** sized for the largest system (14.1). Matrices and vectors are per
system.

**Example.** About 9 systems (u, v, w, p, T, k, eps, plus scalars) at about 280 B per cell each in
matrices and vectors, plus about 1 KB per cell of mesh, gradient side table and scratch, gives
about 3.5 KB per cell before AMG. A 40 GB GPU therefore holds roughly 9M cells after the pressure
AMG hierarchy (2-4x the pressure matrix). This is consistent with the commercial codes' 8-30M cells per GPU on 40-80
GB devices.

---

## 4. Linear system layer

### 4.1 Owned-row COO pattern

Each Kokkos FV linear system with one FV variable gets a fixed COO pattern built once at setup:

```
for each owned cell c (in CellID order) where the variable is active:
    emit (row_global(c), col_global(c))                -> diag_slot[c]
    for side s in 0..n_sides(c)-1:
        f = faces(s, c); n = other cell of f
        if f is an internal face for this variable:
            emit (row_global(c), col_global(n))        -> slot(s, c)
        else:
            slot(s, c) = -1
```

- Rows are owned only. Columns may be ghost-cell global dofs.
- The COO list is sorted by row and has no duplicates, so `MatSetPreallocationCOO` does a trivial
  reduction and the `MatSetValuesCOO` MPI phase is empty. That a pre-sorted list is handled cheaply
  inside PETSc still needs to be verified (Q9).
- **Multiple variables per system** (allowed by the framework, not used by Navier-Stokes) extend
  the emission to a loop over variables, with rows grouped per variable. Phase 0-2 assert one
  variable per Kokkos FV system.

**Relation to the existing Kokkos `Matrix`.** The existing class builds full CSR including ghost
rows and uses `Matrix::find`. FV systems use the new `FVSystemAssembly` class instead (4.2). FE
Kokkos systems are unchanged. Mixed FE/FV systems keep the existing path in phase 0; see Q3.

### 4.2 The `FVSystemAssembly` object

```cpp
class FVSystemAssembly              // one per Kokkos FV linear system; host object
{
  // pattern (built once per mesh change)
  Array<PetscInt>   coo_i, coo_j;   // host only; freed after MatSetPreallocationCOO
  PetscCount        nnz;

  // per-tag device arrays (normally one matrix tag and one vector tag)
  Array<Real>       values;         // [nnz], device; passed to MatSetValuesCOO
  Array<Real>       rhs;            // [n_owned_rows], device; aliases the PETSc RHS local array

  // shared scratch (owned by the problem, sized for the largest system)
  Array2D<Real>     hf_diag, hf_off, hf_rhs;   // (side, c)
  Array<Real>       cell_diag, cell_rhs;       // [n_owned]

  // assembly mode
  enum class Mode { Atomic, Gather } mode;     // R3; Atomic is the default (4.4)
};
```

**As implemented (2026-10-15).** `Moose::Kokkos::FVSystemAssembly`, owned by `LinearSystem`:
- The value array of each matrix tag is laid out as the owned rows of the Kokkos system's sparsity
  pattern, which lists owned rows first, so the slots from `System::getMatrixSlot()` apply
  unchanged. `MatSetPreallocationCOO` gets the owned rows only.
- Right-hand sides are added straight into the local array of each tagged PETSc vector, which
  `LinearSystem` zeroes before assembly; host vectors on a GPU are staged on device.
- Before each dispatch an object gets the device arrays of its tags (up to four of each). The
  Kokkos system no longer activates residual or matrix tags for linear FV, and only provides the
  solution vectors.
- The scratch is per flux kernel (4.4 gather mode), not shared per problem.

### 4.3 Assembly sequence for one linear system

```
A0  zero scratch                          one fused kernel over half-faces and cells
A1  BC relation passes                    only BCs whose relations are not "static" (5.3)
A2  face passes                           one launch per flux kernel over its range list
A3  elemental passes                      one launch per elemental kernel over owned cells
A4  gather (+ fused relaxation)           one launch over owned cells; writes values[] and rhs[]
A5  fence (single)                        before handing device pointers to PETSc
A6  MatSetValuesCOO(values, ADD/INSERT)   device pointer; INSERT since gather writes full rows
```

**Fences (R1).** There are no fences between A0 and A4. All launches go to the default execution
space instance in stream order. A5 is the only fence, and it is needed only if PETSc uses a
different stream from Kokkos.
- PETSc's Kokkos backend can share Kokkos's default execution space instance. If it does, A5
  becomes a no-op. To verify (Q8).
- This requires a framework change: a `Dispatcher` mode with no trailing fence (5.6).

**Time sections.** Use `KokkosFVZero`, `KokkosFVBC`, `KokkosFVFace`, `KokkosFVCell`,
`KokkosFVGather`, `KokkosFVUpload`. Enable the Kokkos Tools NVTX connector for timelines (R12).

### 4.4 Face pass

For a flux kernel `K` on variable `u`, one thread per face in `K`'s range list:

```cpp
template <typename Derived>
KOKKOS_FUNCTION void FluxKernel::operator()(FacePass, ThreadID t, const Derived & k) const
{
  const auto [f, range] = _ranges.locate(t);           // binary search over prefixes
  const FaceDatum face(f, range.orientation, mesh());
  if (range.internal_for_var)
  {
    const auto c = k.template computeInternalFaceContribution<Derived>(face);
    // row E (always owned in INTERIOR/INTERFACE; may be ghost only in PROCESSOR ranges)
    if (face.elemOwned())
    {
      hf_diag(hE) += c.elem;  hf_off(hE) += c.neighbor;  hf_rhs(hE) += c.rhs_elem;
    }
    if (face.neighborOwned())
    {
      hf_diag(hN) += c.neighbor_row_diag();   // -c.neighbor by default (antisymmetric)
      hf_off(hN)  += c.neighbor_row_off();    // -c.elem by default
      hf_rhs(hN)  += c.rhs_neighbor;
    }
  }
  else
  {
    const auto bc = _bc_relations.at(range, f);        // aligned 1:1 with the range (5.3)
    const auto c = k.template computeBoundaryFaceContribution<Derived>(face, bc);
    hf_diag(hLocal) += c.diag;  hf_rhs(hLocal) += c.rhs;
  }
}
```

- **No races.** Every half-face belongs to exactly one face, and kernels run sequentially.
- **The elem and neighbor rows of the 2x2 face block.** Under the CPU convention, row N is
  `(-a_EE, -a_EN)`: its diagonal is `-a_EN` and its off-diagonal is `-a_EE`. The default helpers
  encode exactly that.
- **Non-antisymmetric kernels** override `neighborRowDiag` and `neighborRowOff`, for example
  kernels that skip a block's row, or where a sign issue such as the drift-flux one is found
  intentional.
- **PROCESSOR ranges** check `elemOwned()` and `neighborOwned()` on each thread. In the other
  range kinds these are compile-time or range-uniform, so the branches do not diverge.

**Atomic mode (R3), the default.** The same hook is called. The thread then `atomic_add`s
directly into `values[diag_slot]`, `values[slot]` and `rhs`, skipping the scratch and the gather.
The GPU campaign made this the default (2026-10-14): on an A100 the atomic face pass with fused
terms (V2F) was 2.1x faster than today's scheme, while the gathers were slower than it (V3
0.8-1.0x, V4 0.45-0.73x). Atomic results are not bitwise reproducible.

**Gather mode, opt-in.** For bitwise-reproducible assembly. Phase 0 implements the V3 form, which
was faster than V4 on the GPU: the face pass of each flux kernel writes per-face scratch (the four
entries of the 2x2 face block and the two RHS values), and a cell pass over `elem_faces` sums each
row in local side order and adds it to the tagged matrix and vector. Each entry then receives
exactly one addition per launch, so the existing tag handling is kept. The mode is selected per
problem with `kokkos_linear_fv_assembly = atomic|gather`.

**Fusion.** The navier_stokes kernels already combine their terms in one object (for example
`LinearWCNSFVMomentumFlux` does advection and diffusion), so one face hook per object gives the
benefit of V2F without a composite dispatcher.

### 4.5 Elemental pass

One thread per owned cell where the variable is active:

```cpp
const CellDatum cell(c, mesh());
const auto e = k.template computeCellContribution<Derived>(cell);   // {diag, rhs}
cell_diag[c] += e.diag;  cell_rhs[c] += e.rhs;
```

Elemental kernels: time derivative, sources, Boussinesq, buoyancy, friction, volumetric heat
transfer, pressure gradient, and the TKE/TKED sources.

### 4.6 Gather with fused relaxation

One thread per owned active cell `c`:

```
D = cell_diag[c];   b = cell_rhs[c];   S = 0
for s in sides(c):
    D += hf_diag(s, c);   b += hf_rhs(s, c)
    if slot(s, c) >= 0:  a = hf_off(s, c); values[slot(s, c)] = a; S += |a|
if relax (lambda < 1):                               // equation relaxation, Juretic
    D_star = max(|D|, S) / lambda
    b += (D_star - D) * x_old[c]                     // x_old = previous iterate, device vector
    D = D_star
values[diag_slot[c]] = D;  rhs[c] = b
if record_row_sums:  rowsum[c] = D + sum of off-diagonals   // for SIMPLEC (8.2)
```

- This reproduces `NS::FV::relaxMatrix` and `relaxRightHandSide` (`SegregatedSolverUtils.C:24-109`)
  without host row access.
- The CPU computes `max(|D|, sum |offdiag|)` over the assembled row. Here every contribution to the
  row has been summed before relaxation, so the result is the same.
- Relaxation is skipped when `lambda = 1`, and for pressure. That matches the CPU, where pressure
  has no equation relaxation.

### 4.7 Determinism

In Gather mode, assembly output depends only on the input data and the fixed summation order
(sides in local order, kernels in warehouse order). It is bitwise identical across runs and
across GPUs of the same architecture. PETSc's SpMV, reductions and AMG are not guaranteed to be
deterministic, so whole-solve reproducibility is not claimed.

---

## 5. Object interfaces

### 5.1 Datums

```cpp
class FaceDatum                       // value type, built per thread
{
  FaceID f; int orient; const FVMeshData & m;
public:
  KOKKOS_FUNCTION CellID elem() const;          KOKKOS_FUNCTION CellID neighbor() const;
  KOKKOS_FUNCTION bool elemOwned() const;       KOKKOS_FUNCTION bool neighborOwned() const;
  KOKKOS_FUNCTION Real  area() const;           KOKKOS_FUNCTION Real3 normal() const; // orient applied
  KOKKOS_FUNCTION Real  gC() const;             KOKKOS_FUNCTION Real3 dCN() const;
  KOKKOS_FUNCTION Real  invDn() const;          KOKKOS_FUNCTION Real3 nocVector() const;
  KOKKOS_FUNCTION Real3 centroid() const;       KOKKOS_FUNCTION Real3 skewVector() const;
  KOKKOS_FUNCTION Real  elemVolume() const;     KOKKOS_FUNCTION Real3 elemCentroid() const;
  KOKKOS_FUNCTION FaceID id() const;
};

class CellDatum { CellID c; ... volume(), centroid(), subdomain(), id() };
```

`FVDatum` and `MeshDatum` stay as they are for FE use. The existing Kokkos linear FV objects
migrate to `FaceDatum` and `CellDatum` (Q3).

### 5.2 Flux kernel API

```cpp
struct InternalFaceContribution
{
  Real elem, neighbor, rhs_elem, rhs_neighbor;           // CPU semantics
  KOKKOS_FUNCTION Real neighbor_row_diag() const { return -neighbor; }
  KOKKOS_FUNCTION Real neighbor_row_off()  const { return -elem; }
};
struct BoundaryFaceContribution { Real diag, rhs; };

class Moose::Kokkos::FVFluxKernel : public LinearSystemContributionObject, public BlockRestrictable
{
public:
  struct FacePass {};
  template <typename Derived>
  KOKKOS_FUNCTION InternalFaceContribution computeInternalFaceContribution(const FaceDatum &) const; // default aborts
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryFaceContribution computeBoundaryFaceContribution(const FaceDatum &,
                                                                           const BoundaryRelations &) const;
  // compile-time capability flags, read by the registration macro (as uses_precompute_hooks is today)
  static constexpr bool needs_boundary_value = false;
  static constexpr bool needs_boundary_normal_gradient = false;
};
```

- Registration reuses `registerKokkosLinearFVKernel` with a new SFINAE branch that registers
  `FacePass` when `computeInternalFaceContribution` is user-defined. This mirrors how
  `has_split_linear_fv_flux_dispatchers` selects the loops today
  (`KokkosDispatcher.h:364-435`).
- **Combined matrix and RHS.** One hook returns both. For RHS-only reassembly, a
  `ComputationMode` template parameter lets the face pass skip the matrix writes, and the compiler
  removes the dead matrix code because the mode is a compile-time constant.

### 5.3 BC API and relation storage

- The existing `BoundaryRelation{coefficient, source}` for value and normal gradient is kept
  (`KokkosLinearFVBoundaryCondition.h`).
- **Change:** relation arrays are indexed by the BC's face ranges. The index is
  `prefix[r] + (f - range.begin)`, a contiguous segment per range. This replaces the
  `_kokkos_element_side_ids` `std::set` lookup and also fixes the suspected one-side-per-element
  limitation of `sideWithBoundaryID`. Implemented 2026-10-15: a BC covers the ranges where its
  variable lives on one side, the cell there is owned, and that side carries exactly one boundary
  ID, one of the BC's. A flux kernel stores, per boundary range, the BC index, the offset of the
  range in the BC's arrays and the side.
- **New flags:**
  - `static constexpr bool includes_material_multiplier`, `provides_complete_flux` and
    `needs_boundary_noc`, mirroring the CPU `LinearFVAdvectionDiffusionBC` flags.
  - `bool isStatic() const`: the relations depend only on geometry and constant parameters. Static
    BCs compute their relations once after `initialSetup` and are skipped in A1 (memoization,
    following Albany).
- **NS BCs needing more than a scalar relation:**
  - velocity symmetry needs the other components' cell values;
  - pressure symmetry and pressure flux need `HbyA_f`.

  These read the extra arrays directly in their relation hooks. A BC that depends on other
  equations' values is **not** static.
- **Fusion (R5).** A variable usually has 3-6 BC objects. Today each is its own A1 launch.
  Because most BCs are static, A1 typically launches 0-2 kernels per assembly. Full fusion across
  BC types would need a composite dispatcher and is deferred.

### 5.4 Elemental kernel API

```cpp
struct CellContribution { Real diag, rhs; };
template <typename Derived>
KOKKOS_FUNCTION CellContribution computeCellContribution(const CellDatum &) const;
```

### 5.5 Interpolation methods

These are trivially copyable structs held by value in the kernel:

```cpp
struct AdvectedWeights { Real w_elem, w_neighbor, rhs_face_value; };
struct Upwind        { KOKKOS_FUNCTION AdvectedWeights operator()(const FaceDatum &, Real pE, Real pN,
                                                                  const Real3 & gE, const Real3 & gN,
                                                                  Real flux) const; };
struct Average       { ... };   struct Harmonic { ... };
struct VanLeerWeight { ... };   struct MinmodWeight { ... };   struct MUSCLDeferred { Real factor; ... };

using AdvectedMethod = Variant<Upwind, Average, VanLeerWeight, MinmodWeight, MUSCLDeferred>;
```

- The kernel stores an `enum` plus the parameters.
- The face pass dispatches through a `switch` that is uniform within a launch (all threads take
  the same case). That avoids device virtual calls without multiplying template instantiations.
- The math is ported from `framework/src/fvinterpolationmethods/*` as pure functions.

### 5.6 Required framework changes to the dispatcher

| Change | Why | Where |
|---|---|---|
| A `parallelForNoFence` path, plus an execution-space instance accessor | R1 | `KokkosDispatcher.h:83-86` (currently always `parallel_for` then `fence`) |
| A launch over a range list (prefix and binary search) | 3.6 | new helper in `KokkosThread.h` or the FV base classes |
| Optional `LaunchBounds` per dispatcher | P9 | `Policy` alias in `KokkosDispatcher.h:21` |
| A report of `sizeof(Object)` per registered dispatcher, in debug or under a flag | R10 | `DispatcherRegistry::addDispatcher` |

### 5.7 Functor size policy (R10)

`Dispatcher` stores `const Object _functor_device` by value (`KokkosDispatcher.h:104`), so the
**whole MOOSE object**, including its `MooseObject` and interface base classes, is the kernel
argument.

- With Kokkos' `__grid_constant__` path (Volta+, CUDA 12.1+), objects up to 32 KB are passed as
  parameters.
- Otherwise objects from 512 B to 32 KB go through constant memory with a host-blocking
  serialization point on every launch.

Policy:

1. In phase 0, measure `sizeof` for representative objects, and check which launch path the
   PETSc-downloaded Kokkos 4.7.04 selects on the target toolchain (Q7).
2. If objects exceed the fast path, the fix is framework-wide and needs agreement upstream. The
   option is to split each Kokkos object into a host object and a small device payload struct
   that the dispatcher copies. That changes the Kokkos-MOOSE functor model, so it is not proposed
   unless the measurement demands it.

---

## 6. Field data

### 6.1 Cell values

- A linear FV variable has one dof per cell. With one variable per system, the local dof of cell
  `c` is a direct lookup in the existing `_local_elem_dof_index(0, c)`.
- `CellValue` accessors return `vector[dof]` for the current state, `old`, `older`, or the
  previous iterate (the `solutionPreviousNewton` equivalent).
- The ghost values come from the PETSc ghost update of the solution (6.3).

### 6.2 Gradients

**Storage (R6).** Each variable whose gradient is requested gets a ghosted **block** vector with
block size `dim`: component `d` of cell `c` is entry `dim * dof(c) + d`.
- One ghost update moves all components.
- Old gradient states (needed only for FVReconstructedPressureGradient's state-1 coupling) are
  additional block vectors.
- Double buffering as on the CPU: compute every requested gradient into the "next" vectors, then
  swap pointers.

**Green-Gauss** is a cell gather with no atomics. One thread per owned cell reads a padded
per-(side, cell) side table stored **cell-fastest** (variant G0C of the mini-app):

```
g = 0
for s in sides(c):                                    // side table entry (s, c) at s * n_cells + c
    n = neighbor(s, c)
    if internal for var:  phi_f = gc(s, c)*phi_c + (1 - gc(s, c))*phi_n
    else:                 phi_f = rel.coefficient * phi_c + rel.source   (value relation)
    g += phi_f * area(s, c) * normal(s, c)             // outward normal of c: no sign flip
g /= volume[c];  if RZ: g_r -= phi_c / r_c
```

**Layout decision: G0C (2026-10-14).** The side table holds, per (side, cell), the neighbor, the
boundary ID, the area, the interpolation weight seen from `c`, the outward normal and the face
centroid, with normals and centroids as `Real3`. Each internal face is stored twice, once per
cell.

| Variant (hex 2.1M, A100, ms per gradient) | native | Morton/owner | RCM/owner | scrambled |
|---|---|---|---|---|
| G0: side table, side-fastest (today's Kokkos mesh order) | 1.27 | 1.47 | 1.15 | 2.25 |
| **G0C: side table, cell-fastest** | **0.56** | **0.55** | **0.56** | **0.76** |
| G1C: face arrays through cell-fastest `elem_faces` | 0.48 | 0.51 | 0.55 | 6.02 |
| G1CS: G1C with separate x, y, z arrays | 0.49 | 0.51 | 0.56 | 7.85 |

- **Why G0C over G1C.** G1C is 10-15% faster on well-numbered meshes but 8x slower on the
  scrambled one, because it reads the face data through an index. G0C stays coalesced whatever
  the numbering, so gradients need no renumbering pass and no check of the mesh ordering.
- **Cost.** About 15% more time than G1C on good numberings, and about 432 against 280 B/cell
  (hex) of input data. The side table is in addition to the face arrays of the assembly (3.4), so
  it counts against cells per GPU (P12).
- **Cell-fastest is required.** With the side-fastest order of today's Kokkos mesh, consecutive
  threads read with a stride of the side count: G0 is 2.3x slower than G0C. The CPU prefers
  the opposite order (G0 about 2x faster than G0C on one thread), so the layout is chosen for the
  GPU.
- **Where the table lives.** The per-side arrays of the Kokkos mesh (`_elem_neighbor`,
  `_side_area`, `_side_normal`, ...) are side-fastest and are reached through `KokkosDatum`,
  which FE objects share. Transposing them changes shared code. The linear FV layer therefore
  builds its own cell-fastest side table at mesh setup, and rebuilds it when the mesh changes.
- Gradients are deterministic: each cell sums its sides in a fixed order.

- The Venkatakrishnan limiter is a second cell pass over the same face list.
- Ghost gradients come from the block-vector ghost update. **They are never computed locally on
  ghost cells**, because a ghost cell's face set is incomplete (P14).

**Q4: are PETSc ghosted device vectors usable?**
- If not, the fallback is a Kokkos halo exchange using the system's existing ghost communication
  lists: one packed device buffer per stage, sent with GPU-aware MPI.
- Write the fallback behind the same interface (`exchangeGhosts(fields...)`) so the choice stays
  internal.

### 6.3 Halo exchange plan per SIMPLE iteration

| Stage | Fields exchanged | Exchanges |
|---|---|---|
| After the momentum solves (u, v, w) | u, v, w solutions | 3 (each solve's libMesh localization; see Q4) |
| Momentum gradients (if NOC, deviatoric terms, or limiters) | grad u, grad v, grad w as block vectors | 3, or 1 with a combined `3*dim` block vector |
| `computeHbyA` | `HbyA`, `Ainv` as one block vector of size `2*dim` | 1 |
| After the pressure solve | p | 1 |
| Pressure gradient | grad p | 1 |
| After the velocity correction | u, v, w | 1 (combined block) or 3 |
| Each additional scalar equation | phi, grad phi | 2 |

This gives about 11-14 exchanges per iteration for the flow equations, against more than 20 with
one vector per field.

**Overlap.** Processor ranges come last within each block, so a later optimization can start the
ghost update, run the INTERIOR face pass, wait, then run the PROCESSOR pass. NekRS measured 10-15%
from the same split. Not in the initial phases.

### 6.4 Face fields

- `mdot[f]`, `HbyA_f[f]` and `rhoAinv_f[f]` (`Real3`) are owned by `KokkosRhieChowMassFlux`.
- `mdot` is restartable through `declareRestartableData<Array<Real>>`. Like the CPU face map, it
  depends on the partition.
- They are never exchanged. Processor faces are computed identically on both ranks (D5).
  - **Requirement:** every input to a processor face, meaning cell values, gradients, properties,
    `HbyA` and `Ainv`, must be the owner's value on both ranks. That holds when every such field
    is ghost-updated from its owner and **never** recomputed on ghosts.

### 6.5 Properties

A Kokkos FV property object fills `Array<Real>(n_owned + n_ghost)` per property per state.

| Kind | Evaluation | When |
|---|---|---|
| Constant | none; a `Scalar` is read in kernels | never |
| `KokkosParsedFunction(t, x)` | cell pass over owned and ghost cells (centroids needed for ghosts) | when `t` changes |
| Expression of cell fields (`mu_eff = mu + mu_t`, `rho * cp`) | cell pass over owned cells, then ghost update; or owned and ghost if all inputs are ghosted | once per SIMPLE iteration, or whenever an input field changes |
| Old states | swap at `timestepSetup` | per time step |

**Face values** are interpolated inside kernels: `gC`-weighted average by default, harmonic where
the CPU uses `FVHarmonicAverage`. This matches the CPU `makeCDFace` evaluation for variable-backed
and smooth properties. For *functor materials whose face evaluation is not an interpolation*, such
as the mixture material, which evaluates its component properties at the face, results can differ
from the CPU at round-off to truncation level. Note this when doing parity testing.

**As implemented (phase 1 step 7, 2026-10-15).** Properties match the CPU exactly rather than to
truncation level, because the CPU semantics turned out to be specific: a functor material is
interpolated from its cell values only on faces evaluated without a side (internal faces through
`makeCDFace()`), and evaluated with the face values of its inputs on faces evaluated from one side
(boundary faces through `singleSidedFaceArg()`); functions are evaluated at face centroids, with a
1e-8 offset on one-sided block-interface faces. `Moose::Kokkos::LinearFVPropertyInterface` mirrors
these rules; expressions of variables use `KokkosParsedFunctorMaterial`. Constant properties are
scalars, function properties are refilled when the time changes, variable-based ones before every
assembly. See `phase1_implementation_plan.md`, step 7.

---

## 7. Navier-Stokes objects

| Kokkos object | Pass | Reads | Writes | Phase |
|---|---|---|---|---|
| `KokkosLinearWCNSFVMomentumFlux` | face | `mdot`, `mu`/`mu_eff` cells, u_i cells, grad u_i, grad u_j (deviatoric), `noc`, advected method | half-face scratch | 2 |
| `KokkosLinearFVMomentumPressure` | cell | grad p (coupling state) | cell scratch (RHS) | 2 |
| `KokkosLinearFVPressureCorrectionDiffusion` | face | `rhoAinv_f`, grad p (NOC/tangential) | half-face scratch | 2 |
| `KokkosLinearFVDivergence` | face (RHS only) | `HbyA_f` | half-face RHS | 2 |
| `KokkosLinearFVTimeDerivative` | cell | rho (current, old), u old states | cell scratch | 3 |
| `KokkosLinearFVEnergyAdvection`, `KokkosLinearFVScalarAdvection`, `KokkosLinearFVDiffusion` | face | `mdot`, properties, gradients | half-face scratch | 3 |
| Boussinesq, buoyancy, friction, volumetric heat transfer, source | cell | properties, other variables | cell scratch | 3 |
| BCs: Dirichlet, outflow, extrapolated, extrapolated pressure, velocity symmetry, pressure symmetry, normal velocity, pressure flux | BC relation | functions, component cells, `HbyA_f`, `rhoAinv_f` | relation arrays | 2-3 |
| Turbulence: TKE/TKED sources, turbulent advection and diffusion, wall-function BC, `mu_t` aux | cell/face | wall data arrays (7.1) | scratch, aux vector | 4 |
| CHT, two-phase, reconstructed pressure gradient | various | | | 4 |

**Shared physics functions** live in a header such as `NSFVKokkosFunctions.h`:

- `pressureDiffusionFace(face, rhoAinv_f, grad_p_E, grad_p_N, noc_flag)` returns
  `{a_E, a_N, b_E}` per unit area. It is used by both the pressure diffusion kernel and the mass
  flux face pass, replacing the CPU's reuse of a mutated kernel object.
- `momentumStressFace(...)` and `advectedWeights(...)`.

### 7.1 Turbulence wall data

The CPU keeps `std::map<const Elem *, ...>` wall maps and allocates vectors per cell. On device:

- `wall_bounded[c]` is a `uint8_t` per cell.
- `wall_face_offset[c]` and `wall_faces[k]` form a CSR list of wall faces for the wall-bounded
  cells only (most cells have none). `wall_dist[k]` holds the matching distances.
- y+ comes from a fixed-iteration device Newton loop. The CPU iteration limit and tolerance become
  compile-time constants, so the loop length is bounded.
- Skipping the turbulent advection and diffusion rows for wall-bounded cells moves into the
  gather: if `wall_bounded[c]`, the row comes only from the TKE/TKED elemental pin. This replaces
  the CPU kernels' direct matrix writes that bypass tags.

### 7.2 Porting inventory

Everything in `modules/navier_stokes` that a GPU-resident linear FV solve needs, in dependency
order. Groups 0 and 1 are a usable laminar SIMPLE/PIMPLE solver; groups 2 and 3 add physics. The
groups map onto the phases of section 14 (group 0 is phases 0-1, group 1 is phase 2, group 2 is
phase 3, group 3 is phase 4).

**Group 0: framework prerequisites.** Every navier_stokes kernel depends on these. Today's Kokkos
linear FV objects lack all of them.

| Piece | Today (Kokkos linear FV on the branch) | Needed for | Design |
|---|---|---|---|
| Cell gradients on device (Green-Gauss, optional limiter) | None; kernels never read gradients | Pressure gradient, non-orthogonal correction (NOC), deviatoric stress, higher-order advection | 6.2 |
| Property arrays (`rho`, `mu`, `cp`, `k`, ...) | Parsed functions only | Every navier_stokes kernel; replaces functor materials and the virtual calls of fluid-property objects | 6.5 |
| Advection with a per-face flux and interpolation methods | Constant velocity, upwind only | Momentum, energy and scalar advection with the Rhie-Chow mass flux | 5.5 |
| Time derivative | Test-only object | PIMPLE and transient SIMPLE | 7 |
| Reading variables of other systems | None | Velocity components in momentum (deviatoric terms, symmetry BC), velocity data in pressure | 3.6 |
| Framework BC set | Functor Dirichlet and Neumann only | Robin, extrapolated, outflow, symmetry | 5.3 |
| Non-orthogonal correction | None (needs gradients) | Diffusion and momentum on skewed meshes | 4.4, phase 1 |

**Group 1: laminar incompressible SIMPLE/PIMPLE.**

| CPU object | Role | What changes on device |
|---|---|---|
| `LinearWCNSFVMomentumFlux` | Momentum advection with `mdot`, viscous terms, NOC, optional deviatoric terms | Face pass reading `mdot`, `mu`, u_i and all velocity gradients (section 7 table) |
| `LinearFVMomentumPressure` | `-dp/dx_i V` source in the momentum RHS | Cell pass reading the coupling pressure gradient |
| `LinearFVPressureCorrectionDiffusion` | Pressure equation with the per-face `rho Ainv` coefficient | Face pass using `pressureDiffusionFace` (section 7) |
| `LinearFVDivergence` | RHS divergence of `HbyA` | Face pass, RHS only |
| `LinearFVExtrapolatedPressureBC`, `LinearFVNormalVelocityFunctorDirichletBC`, `LinearFVVelocitySymmetryBC`, `LinearFVPressureSymmetryBC`, `LinearFVPressureFluxBC` | Flow BCs | BC relation arrays (5.3); symmetry and pressure-flux BCs read other components, `HbyA_f` and `rhoAinv_f` |
| `RhieChowMassFlux` | Face mass flux, `HbyA`, `Ainv`, velocity correction | `KokkosRhieChowMassFlux` (8.1). The face maps keyed by `FaceInfo::id()` become face arrays. `computeHbyA`, `populateCouplingFunctors` and `computeFaceMassFlux` become device passes (8.2 steps 8, 9, 13); today they are serial host loops, and the last two call kernel and BC hooks once per face. `updateCellVelocity` is already vector operations (step 16). |
| `SegregatedSolverUtils` | Equation and field relaxation, normalization factor, pressure pin | The row-sum loops (`get_row`) in equation relaxation and SIMPLEC become the fused relaxation and recorded row sums of the gather (4.6). The normalization factor (`MatMult`), field relaxation and pressure pin are already PETSc operations on device vectors. |

**Group 2: energy and scalar transport.**

- Kernels: `LinearFVEnergyAdvection`, `LinearFVScalarAdvection`, `LinearFVVolumetricHeatTransfer`.
- `LinearFVEnthalpyFunctorMaterial`. Its lambdas become property expressions (6.5).
- `LinearFVConvectiveHeatTransferBC`.
- Conjugate heat transfer: `LinearFVRobinCHTBC`, `LinearFVDirichletCHTBC`, and the `CHTHandler`
  fixed-point loop over interface faces, with its face maps.
- The aux recovery of `T` from `h`, which runs every iteration (8.6).

**Group 3: turbulence, multiphase, porous and buoyancy.**

- **k-epsilon:** `LinearFVTurbulentAdvection`, `LinearFVTurbulentDiffusion`,
  `LinearFVTKESourceSink`, `LinearFVTKEDSourceSink`, `LinearFVTurbulentViscosityWallFunctionBC`,
  and the `mu_t` functor materials (`INSFVkEpsilonViscosityFunctorMaterial`,
  `MixingLengthTurbulentViscosityFunctorMaterial`). These are the hardest objects: per-element
  `std::map` wall data, per-cell allocations, iterative y+ solves, and matrix writes that bypass
  tags. See 7.1.
- **Two-phase mixture:** `LinearWCNSFV2PMomentumDriftFlux`, `NSFVMixtureFunctorMaterial`,
  `LinearFVPhaseChangeSource`.
- **Other sources:** `LinearFVMomentumFriction` (`LinearFrictionFactorFunctorMaterial`),
  `LinearFVMomentumBuoyancy`, `LinearFVMomentumBoussinesq`, `LinearFVRZViscousSource`.
- `FVReconstructedPressureGradient`, with its dense least-squares solve per cell.

**Not ported.**

- **The executioners** (`SIMPLE`, `PIMPLE`, `LinearAssemblySegregatedSolve`, `SIMPLESolve`,
  `PIMPLESolve`). Their control flow stays on the host, and every step they call runs on the
  device (8.6). Q5 asks whether the Kokkos-aware solve is a mode of
  `LinearAssemblySegregatedSolve` or a separate class.
- **The `WCNSLinearFV*Physics` actions.** They only need to create the Kokkos objects.
- **Fluid-property user objects.** These are a non-goal (1.3).

---

## 8. Rhie-Chow and SIMPLE on device

### 8.1 `KokkosRhieChowMassFlux`

A Kokkos general user object, `execute_on = NONE`, driven by the executioner. It holds the face
fields (6.4) and the cell block vector `[HbyA | Ainv]`.

### 8.2 One SIMPLE iteration, step by step

The table lists data movement and synchronization for each step. `F` is a fence, `X` a halo
exchange, `R` a scalar device-to-host read.

| # | Step | Implementation | F/X/R |
|---|---|---|---|
| 1 | Snapshot grad p for coupling | device copy into the coupling-state vector | |
| 2 | Momentum gradients | cell gathers for u, v, w into one block vector | X1 |
| 3 | Assemble the momentum matrix | A0-A4 for the u system, including the momentum-pressure source and fused relaxation | F1 |
| 4 | Share across components (R7) | if the velocity BCs are component-independent (8.3), assemble only the RHS for v and w (`ComputationMode::RHS`) and reuse u's matrix values | |
| 5 | Normalization factors | 2 `MatMult` + vector operations per component; batch the 3 norms into one reduction where PETSc allows | R (x3 or x1) |
| 6 | Solve u, v, w | one KSP, preconditioner set up once, three solves; tolerances looser on non-final outer iterations | R (KSP) |
| 7 | Ghost update u, v, w | libMesh localization of each solution | X (3) |
| 8 | `computeHbyA` | per component: `MatGetDiagonal` -> `A`; `MatMult(M, u) - A.*u - b - gradp_snapshot.*V` -> `H`; `HbyA = H./A`; `Ainv = V./A`; SIMPLEC uses the recorded row sums (4.6); then pack into the block vector | X1 |
| 9 | Face coupling fields | face pass over flow ranges: `HbyA_f`, `rhoAinv_f`; boundary ranges use velocity BC relations, with the CPU Dirichlet special case | |
| 10 | Assemble pressure | A0-A4 for p: pressure diffusion + divergence, no relaxation, pin via diagonal and RHS edit in the gather | F1 |
| 11 | Solve p | CG/AMG with the freeze policy (8.4) | R (KSP) |
| 12 | Ghost update p | localization | X1 |
| 13 | Mass flux | face pass: `mdot = -HbyA_f + pressureDiffusionFace(p)` per unit area times area | |
| 14 | Relax pressure | `VecAXPBY` against the previous iterate | |
| 15 | grad p | cell gather into the block vector | X1 |
| 16 | Velocity correction | cell pass: `u_i = -HbyA_i - Ainv_i * dp/dx_i` | X1 (u, v, w block) |
| 17 | Energy, scalars, turbulence | each: gradients (X), A0-A4 (F), norm (R), solve (R), localization (X), optional field relaxation and limiter as cell passes | per equation |
| 18 | Per-iteration `execute(EXEC_NONLINEAR)` | Kokkos aux kernels and user objects only (R9) | |
| 19 | Convergence | compare the per-equation normalized initial residuals already read in steps 5, 6, 11 and 17 | |

Per flow iteration this gives about 3 fences, about 8-10 halo exchanges, and about 10 scalar
reads outside the KSPs.

### 8.3 Sharing the momentum matrix (R7)

- Valid when, for every velocity BC on every boundary, the matrix part of the relation (value and
  gradient coefficients) is identical across components. Dirichlet, no-slip, outflow and
  extrapolated BCs qualify. Symmetry (`1 - n_i^2`) and normal-velocity BCs do not.
- **Detection:** at setup, each velocity BC reports `componentIndependentMatrix()`.
- **When valid:**
  - assemble the full system for u;
  - assemble RHS only for v and w;
  - `computeHbyA` uses the one matrix: `Ainv` is identical across components, and `H` is computed
    per component with that component's solution and RHS.
- **When invalid:** assemble three matrices as the CPU does.
- MOOSE keeps three `LinearSystem`s either way, which preserves output, restart and
  postprocessing. The shared matrix is a `Mat` reference held by the executioner. PETSc solves the
  three right-hand sides sequentially with one KSP and its preconditioner (Nalu-Wind's
  `HypreUVWSolver` pattern).

### 8.4 Pressure solver and reuse (R8)

- **Default:** CG + BoomerAMG with hypre ≥ 2.31 and the PETSc GPU defaults (PMIS, ext+i,
  l1-Jacobi, keep-transpose), plus aggressive coarsening on 2 levels, `P_max 3`, `truncfactor 0.1`,
  and a strong threshold tuned over 0.25-0.6.
  - Enable a hypre device memory pool.
  - Option sets are in `gpu_fv_design_research.md` section 5.
- **Freeze policy.** Keep the preconditioner while the iteration count stays within
  `1.5 x its_after_setup`. Rebuild otherwise, or every `pressure_pc_rebuild_frequency` outer
  iterations (default 10, to tune), or at the first outer iteration of a time step.
  - This extends the existing `pressure_pc_recompute_frequency` parameter of the CPU executioner.
- **Alternative to benchmark:** GAMG with `-pc_gamg_reuse_interpolation true` (coefficient-only
  Galerkin resetup on device).
- **Momentum and scalars:** BiCGStab or GMRES(30) + Jacobi, or hypre single-level l1-Jacobi.
  Relative tolerance 1e-2 to 1e-3 on non-final iterations.

### 8.5 Norms and convergence

- The CPU convergence test uses each system's KSP *initial* residual divided by its normalization
  factor. Keep this exactly, so convergence histories are comparable with the CPU.
- That makes 2 scalar reads per equation per iteration. That is acceptable: it is roughly 10 reads
  per iteration, against hundreds of kernels.

### 8.6 Host role and device residency

**Rule.** Every per-cell, per-face and per-nonzero computation of a SIMPLE or PIMPLE iteration,
and all the data it touches, is on the device. Only the control flow stays on the host: which
step runs next, how many corrector passes, when to stop.

**The control flow cannot move to the device:**

1. **The linear solves are `KSPSolve` calls.** With PETSc's Kokkos backend, every vector and
   matrix operation of the Krylov method runs on the device. The Krylov loop and its convergence
   test are host code, and no device kernel can call PETSc.
2. **The steps need global synchronization.** Each Krylov iteration, and each SIMPLE step that
   depends on the previous one, needs a reduction over all cells and, with several GPUs, over MPI.
   A kernel cannot wait for the whole device or call MPI, so the host sequences the launches.
3. **Convergence decisions need a scalar on the host.** Stopping SIMPLE, or the corrector count of
   PIMPLE, depends on residual norms, so each check costs one small device-to-host read.

GPU CFD codes (PETSc's and AmgX's GPU paths, SPUMA, NekRS) are structured the same way. What
matters is that the host only launches kernels and reads scalars.

**What the host does per iteration:**

- launch the assembly, relaxation, gradient, Rhie-Chow and correction kernels;
- call `KSPSolve` on device-resident matrices and vectors;
- read the scalars in 8.2 and 8.5: normalization factors and KSP initial residuals, about 10 per
  flow iteration.

**Host work today that must move to the device.** Each of these is a serial host loop, or a
mesh-sized host-device copy, every iteration.

| Today | On device |
|---|---|
| `RhieChowMassFlux::computeHbyA`, `populateCouplingFunctors`, `computeFaceMassFlux` (serial host loops, face hash maps) | Device passes, 8.2 steps 8, 9, 13 |
| `SegregatedSolverUtils` row sums through `get_row` (equation relaxation, SIMPLEC) | Fused into the gather (4.6) |
| Gradient computation and the pressure-gradient snapshot | Block-vector gradients (6.2) and a device copy (8.2 step 1) |
| Per-iteration aux kernels and user objects (`T` from `h`, `mu_t`, wall functions) | Kokkos aux kernels and user objects only (section 10) |
| Kokkos assembly copying solution vectors in and the matrix out through host memory | Device PETSc types, so these are device-to-device copies (implemented on the branch for GPU builds) |

**Reducing the remaining host involvement:**

- **Fewer launches:** fuse the face passes, and capture fixed launch sequences between MPI points
  with `Kokkos::Graph` (13.4).
- **Fewer scalar reads:** write the residuals of all equations into one device array and read it
  once per iteration. Optionally check SIMPLE convergence only every k iterations, if parity
  tests show the outer iteration count does not change materially.
- **Fewer reductions inside the solves:** Krylov methods with fewer global reductions (BiCGStab,
  pipelined variants). The GPU campaign measured BiCGStab with Jacobi 1.7x faster than GMRES(30),
  and GMRES did not converge at 2.1M cells within 500 iterations.
- **PIMPLE:** the number of correctors is fixed (`num_piso_iterations`), so the correctors add no
  convergence reads. The time-step loop adds only the scalars that time stepping needs.

**Verification.** No mesh-sized host traffic in the iteration loop (target in 1.2). Check with
Nsight Systems or rocprof, and on CPU builds with the event counter
(`bench/tools/kokkos_event_counter.cpp`), which reports allocations, fences and launches by
calling function.

---

## 9. Parallel design

| Topic | Decision |
|---|---|
| Ranks per GPU | 1. Document it; warn when `n_ranks_per_node > n_gpus_per_node` and Kokkos FV systems exist. |
| Ghosting | One FV ghost layer, as today. No second layer: ghost gradients are exchanged, not recomputed. |
| Processor faces | Computed on both ranks, each writing its own rows. Bitwise-identical face fluxes require owner-consistent inputs (6.4). Debug builds can check this with a periodic cross-rank comparison of processor-face `mdot`. |
| Partitioning | Standard cell-balanced partitioners. For CHT, consider weighted partitioning later (P15). |
| GPU-aware MPI | Configurable and measured per machine. Device buffers only, never managed memory. |
| Off-process assembly | None (owned rows only). The `MatSetValuesCOO` MPI phase is empty. |

---

## 10. CPU object policy (R9)

At `initialSetup` of a Kokkos FV system:

- **Error** if any non-Kokkos linear FV kernel or BC is attached to the system. Kokkos flux
  kernels already reject non-Kokkos BCs (`KokkosLinearFVFluxKernel.K:187-213`).
- **Error** if a Kokkos kernel's property or field input is provided by a non-Kokkos object
  (functor material, CPU aux kernel). Without this check, the input would silently force host
  evaluation and copies every assembly.
- **Warn** if any non-Kokkos object executes on `NONLINEAR` or `LINEAR`, because
  `FEProblemBase::computeLinearSystemTags` runs those on every assembly.
- **Warn** if the PETSc matrix or vector types of the system are host types.

---

## 11. Instrumentation (R12)

- **Time sections** for each step in 8.2 and each sub-step in 4.3.
- **Kokkos Tools.** Use `space-time-stack` for kernel counts and fences, and the NVTX connector
  for Nsight Systems timelines. Document the environment variables in the docs page.
- **Counters printed with `-log_view`-style output:**
  - kernel launches per SIMPLE iteration;
  - fences per SIMPLE iteration;
  - halo exchanges per iteration;
  - bytes of halo traffic;
  - KSP iterations per equation;
  - preconditioner rebuilds.
- **Rule from the research:** profile assembly against the solves before tuning either. SPUMA found
  assembly small and the solves dominant after porting.

---

## 12. Verification plan

### 12.1 Unit-level (device and host)

| Test | Checks |
|---|---|
| Face construction | Each internal face appears once. The two half-faces of an internal face agree. The sum over a closed cell of `sign * normal * area` is about 0. `gC`, `dCN` and `area` match `FaceInfo` exactly. Range tables cover every face exactly once. |
| Renumbering | A permutation (bijection). Block ranges are preserved. FE Kokkos tests give unchanged results within tolerance. |
| COO pattern | Equals the libMesh sparsity restricted to owned rows. Slots are unique. Sorted by row and column. |
| Gather | Atomic and gather modes produce the same matrix to round-off. Running gather twice is bitwise identical. |
| Relaxation | Matches `NS::FV::relaxMatrix` and `relaxRightHandSide` on CPU-assembled matrices to round-off. |
| Green-Gauss | Linear fields have exact gradients on orthogonal meshes and match the CPU on skewed meshes. |

### 12.2 Parity against the CPU

For each ported kernel and BC, assemble the same system with the CPU and Kokkos objects and
compare the matrices and right-hand sides entrywise. Use relative tolerances around 1e-12, and
looser ones where the face-value evaluation of properties differs (6.5).

### 12.3 Regression

- Phase 0: `test/tests/kokkos/linearfvkernels/*` and `test/tests/kokkos/fe-fv-interplay`.
- Phase 2: `modules/navier_stokes/test/tests/finite_volume/ins/{lid-driven,channel-flow,mms}/linear-segregated`.
  - Iteration histories should match the CPU to about 1e-8 relative per iteration.
  - The MMS input must keep its convergence order.
- Phase 3-4: `block_restriction/linear-segregated`, the turbulence inputs, and the transient
  variants.
- Run every test with `-p 1`, `-p 3` (processor faces) and `--recover` (restart of `mdot`).
  Following `AGENTS.md`, do not create separate specs for these modes.

### 12.4 Performance benchmarks

- **Cases.** A 3D lid-driven cavity and a channel, at 1M, 4M and 16M cells per GPU, with 1, 2, 4
  and 8 GPUs (strong and weak scaling).
- **Compare against** the CPU at equal node count, and against a full CPU node.
- **Metrics:**
  - time per SIMPLE iteration, split by step (8.2);
  - the counters in section 11;
  - achieved bandwidth of the face pass and gather against STREAM. miniAero's best reached 34-36%.

---

## 13. Extensibility notes

### 13.1 Multiple variables per system

The COO builder and gather generalize by emitting rows per variable. FV kernels would then need a
variable-local half-face scratch, so the scratch grows by the variable count. Deferred.

### 13.2 Writing directly into PETSc's values

Instead of `values[]` plus `MatSetValuesCOO`, write into the PETSc AIJKOKKOS CSR value array
(`MatSeqAIJGetCSRAndMemType` on the diagonal block, plus the off-diagonal block for ghost
columns). That saves one array and one pass per system, but couples us to PETSc internals. Do it
only if the upload shows up in profiles.

### 13.3 Block-coupled pressure-velocity

The face pass produces per-face 2x2 blocks per variable pair. A coupled solver would emit 4x4
blocks into a block COO pattern built by the same machinery. Nothing in sections 3-6 prevents
this.

### 13.4 Kernel fusion and graphs

After phase 2 profiling:

1. Fuse the face passes of all flux kernels of one variable through a composite registered at the
   system level.
2. Capture fixed sequences between MPI points with `Kokkos::Graph`. The research shows this pays
   only for sequences replayed at least about 40 times with small kernels, and that MPI and PETSc
   calls cannot be inside a graph.

---

## 14. Phased plan

| Phase | Deliverables | Exit criteria |
|---|---|---|
| **0 - Infrastructure** | Face SoA, cell SoA, renumbering (off by default), range tables, half-face mapping, COO builder, `FVSystemAssembly` with gather and atomic modes, no-fence dispatch, range-list launches, `FaceDatum`/`CellDatum`, time sections and counters, functor-size report. Migrate the existing Kokkos linear FV diffusion, advection, source and functor Dirichlet/Neumann objects. | 12.1 passes. Existing Kokkos linear FV tests give unchanged results. One fence per assembly. Functor sizes measured and Q7 answered. Q4 and Q8 answered by prototype. |
| **1 - FV numerics** | Cell value accessors, block-vector gradients (Green-Gauss and limiter), NOC, interpolation-method structs, constant and parsed property arrays, static BC memoization. | Kokkos diffusion with NOC and limited advection match the CPU to round-off on skewed meshes (12.2). |
| **2 - SIMPLE** | `KokkosRhieChowMassFlux`, momentum flux, momentum pressure, pressure diffusion and divergence, flow BCs, Kokkos-aware `LinearAssemblySegregatedSolve` with fused relaxation, momentum matrix sharing, pressure reuse policy, CPU object policy checks, default device PETSc types. | 12.3 phase-2 tests match the CPU. 12.4 benchmarks run. No mesh-sized host traffic in the loop (verified in Nsight). |
| **3 - Physics breadth** | PIMPLE, time derivative, energy, scalars, Boussinesq and buoyancy, pressure-flux BC, property expressions. | Corresponding CPU tests. |
| **4 - Advanced** | k-epsilon (wall data arrays, device wall functions), CHT, two-phase, reconstructed pressure gradient. | Corresponding CPU tests. Verify the drift-flux sign issue first. |
| **5 - Optimization** | Interior/processor overlap, kernel fusion, graph capture, writing directly into PETSc values, renumbering on by default if it helps. | Benchmark-driven. |

---

## 15. Risks and open questions

### 15.1 Top risks

| Risk | Mitigation in this design | Residual |
|---|---|---|
| Launch and fence latency dominates below about 2M cells per GPU (P1) | One fence per assembly; static BC memoization; combined matrix and RHS hooks; later fusion and graphs | The SIMPLE structure inherently has many steps |
| Linear solves dominate after porting (P4/P5) | Momentum sharing, AMG freezing, memory pool, loose non-final tolerances | AMG setup on GPU; outer iteration count (P17) |
| Halo latency (P2) | Block vectors, about 11-14 exchanges per iteration | Depends on Q4 and on the machine's MPI stack |
| Functor size forces the constant-memory launch path (P10) | Measure in phase 0 | A framework-wide change if it is bad |
| Memory per GPU limits cells per GPU (P12) | Shared scratch, owned rows only, optional face fields | PETSc's and AMG's own storage |
| Parity differences in property face values (6.5) | Documented; tolerance in parity tests | Users comparing CPU and GPU runs |
| Upstream coordination on replacing the existing Kokkos FV path (Q3) | Phase 0 migration keeps the existing tests passing | Needs agreement from the authors of #32492 |

### 15.2 Open questions

Q1-Q6 are carried over from the overview. Q1 (flat ranges) and Q2 (build from `FaceInfo`) are
assumed resolved as recommended there. Q7-Q10 were added with this design, and Q11-Q19 on
2026-10-14 for the GPU campaign.

**Measurement questions** are answered by `bench/run_gpu_campaign.sh`; the section names are in
the last column. CPU results to date are in `findings_2026-10-14.md` and
`bench/results/GRADIENT_RESULTS_2026-10-14_openmp_m4max.md`.

| Q | Question | Status or CPU evidence | Campaign section |
|---|---|---|---|
| Q4 | Are ghost updates of libMesh's ghosted vectors with Kokkos vector types done on the device? Decides between PETSc vectors and the custom exchange (6.2). | Open | `q4` (`-log_view` copy counts on 2-4 ranks) |
| Q7 | Size of Kokkos linear FV functors, and the CUDA launch path (5.7) | **Answered.** Dispatched objects are 4.8-7.5 KB (advection 4824 B, BCs 5408 B, diffusion 6440 B, source 7536 B). Kokkos CUDA passes functors under 512 B as kernel arguments and stages 512 B-32 KB through constant memory, so every linear FV launch takes the constant-memory path. Its cost is Q14. | (Q14) |
| Q8 | Does PETSc's Kokkos backend run on Kokkos's default execution space instance? If so, ordering does not need the fence before `MatSetValuesCOO`. | Open on GPUs. Trivially yes on OpenMP. | `q8` (`tools/petsc_kokkos_probe.cpp`) |
| Q9 | Is `MatSetValuesCOO` on AIJKOKKOS cheap for owned-row, sorted, duplicate-free COO lists? | Open | `overhead` (`kokkos_upload` section), `q4` (`MatSetValuesCOO` row) |
| Q11 | Gradient layout: padded per-(side, cell) side data (G0, G0C) or the face array through `elem_faces` (G1, G1C, G1CS); gather or atomics (G2, G3) | **Decided: G0C, a gather** (6.2). A100, hex 2.1M: G0C 0.55-0.56 ms on good numberings and 0.76 ms scrambled; G1C 0.48-0.55 ms but 6.02 ms scrambled; G0 1.15-2.25 ms. The atomic variants G2 and G3 are slower and not reproducible. CPU: gathers beat atomics by 1.5-1.8x. | `miniapp` (gradient sweep) |
| Q12 | Face normals and centroids as `R3` (AoS) or separate x, y, z arrays (SoA); per-cell tables side- or cell-fastest | **Answered: `Real3`, cell-fastest.** On the A100, cell-fastest makes the side table 2.3x faster (G0C against G0); separate x, y, z arrays give nothing on good numberings and are slower scrambled (G1CS 7.85 against G1C 6.02 ms). The CPU prefers side-fastest. | `miniapp` (G1 vs G1C vs G1CS, G0 vs G0C) |
| Q13 | Fusing the face passes of all flux kernels of a variable (13.4) | CPU: the fused pass (V2F) is 1.5-1.65x faster than V2 at 1.5-2.1M cells (1 and 8 threads, native and RCM ordering), with identical results to 7e-16 | `miniapp` (layout sweep) |
| Q14 | Launch cost, and the cost of the constant-memory launch path for 6 KB functors | CPU (OpenMP): about 9 us per launch; the 6 KB functor adds about 1 us | `miniapp` (`--launch-bench`) |
| Q15 | Cost of the parsed-function interpreter on device (divergence, registers) | CPU: replacing the parsed diffusion coefficient with a constant cuts assembly by about 30% at 4k cells | `coefficients` |
| Q16 | Linear solvers on device: CG with Jacobi, BoomerAMG or GAMG with reuse for Poisson; GMRES, BiCGStab or pipelined BiCGStab for advection-diffusion (8.4) | Earlier GPU run: BiCGStab with Jacobi 1.7x faster than GMRES(30) | `solvers` |
| Q17 | Ranks per GPU and multi-GPU scaling, including halo cost (section 9) | Open | `ranks` |
| Q18 | Cell numbering and face order on the GPU (native, RCM or Morton with owner-sorted faces, scrambled) | A100: RCM and Morton beat native by at most 3% in assembly; scrambled numbering costs assembly up to 3.6x, G1C 12x and G0C 1.4x. Keep native numbering; with G0C, gradients do not need renumbering. CPU: scrambled numbering costs G1 7-12x. | `miniapp` (both sweeps) |
| Q19 | Per-assembly host overhead after the fence and allocation changes, and the default device types | CPU: allocations per assembly 26 to 3, fences about 220 to 54. Earlier GPU run: about 2 ms per assembly before. | `overhead`, `types` |

**Design decisions** (not measurements):

- **Q3.** Replace or coexist with the existing Kokkos linear FV `(side, elem)` path. Q11 and Q12
  inform it.
- **Q5.** Should the Kokkos-aware solve be a mode of `LinearAssemblySegregatedSolve` or a separate
  class?
- **Q6.** Target hardware for the benchmarks: NVIDIA only, or AMD as well. The campaign runs on
  either.
- **Q10. Resolved:** an error by default. A linear system that mixes Kokkos and non-Kokkos linear
  FV objects errors unless `Problem/allow_mixed_kokkos_linear_fv = true`; implemented on the
  branch.
