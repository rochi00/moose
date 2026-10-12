# Kokkos Navier-Stokes Linear FV: Proposed Design

Design proposal for a GPU port of the Navier-Stokes linear FV solver (SIMPLE/PIMPLE with
Rhie-Chow) on top of Kokkos-MOOSE. It builds on:

- `kokkos_moose_design_choices.md`: Kokkos-MOOSE conventions.
- `ns_linearfv_design.md`: the existing CPU Navier-Stokes linear FV design.

References are to MOOSE `next` at `410cc90c7e`.

## Goals

In priority order:

1. **GPU scalability and speed.** This is the main goal.
2. **Alignment with the existing Kokkos-MOOSE design.**
3. **A vectorized (structure-of-arrays) container for each MOOSE mesh block.**
4. **Faces that store indices into an element array, not pointers.** The face storage is itself a
   flat array.
5. **For every element, the indices of the faces that surround it.**

---

## 0. Summary of key decisions

| # | Decision | Goal |
|---|---|---|
| D1 | Add a **face index space** to `Moose::Kokkos::Mesh`. Each mesh face is stored once, in structure-of-arrays (SoA) form, indexed by a dense `FaceID`. | 1, 4 |
| D2 | A face stores `elem[f]` and `neighbor[f]` as **contiguous element IDs** into the element SoA arrays. There are no pointers. | 4 |
| D3 | **Element-to-face connectivity** `elem_faces(side, e) -> FaceID`, padded to the maximum side count, laid out so each element's faces are contiguous. | 5 |
| D4 | Order faces **by block, then by kind** (interior, processor-boundary, block-interface, domain boundary per boundary ID). Each group is a contiguous range, so a block's "container" is a slice of the flat arrays. | 3, 1 |
| D5 | Orient faces with the **CPU `FaceInfo` ownership rule**, which uses global IDs. Both ranks then compute identical fluxes on processor-boundary faces, with no face communication. | 1 |
| D6 | Each rank writes **only its own rows**. Processor-boundary faces are computed redundantly on both sides. Off-process matrix/RHS assembly is never needed. | 1 |
| D7 | Assembly in two steps: a **face pass** computes each face's coefficients once; a **cell gather pass** forms diagonals and right-hand sides. Off-diagonal entries go into precomputed matrix slots. No atomics, no searches, deterministic. | 1 |
| D8 | Keep Kokkos-MOOSE idioms: templated `KOKKOS_FUNCTION` hooks, the dispatcher registry, `registerKokkos*` macros, holders, value binding, `Kokkos<CPUName>` objects in the normal input blocks. | 2 |
| D9 | Store face fields (mass flux, interpolated `HbyA`, `rho * Ainv`) as **face-indexed device arrays**. Store cell fields (gradients, `HbyA`, `Ainv`) as PETSc-backed vectors aliased on device. | 1, 2 |
| D10 | Replace polymorphic functors with **precomputed cell and face property arrays**, filled by Kokkos objects before assembly. | 1 |
| D11 | Keep the SIMPLE/PIMPLE linear algebra in **PETSc device objects**. Fuse equation relaxation into the assembly gather pass. | 1 |

---

## 1. Alignment with Kokkos-MOOSE (goal 2)

### Kept as-is

- `.K` sources, the `MOOSE_KOKKOS_SCOPE` / `MOOSE_KOKKOS_ENABLED` split, and no relocatable
  device code (RDC). So device functions live in headers.
- Objects named `Kokkos<CPUName>` (`KokkosLinearWCNSFVMomentumFlux`, `KokkosRhieChowMassFlux`,
  ...) in `Moose::Kokkos`, placed in the normal `[LinearFVKernels]`, `[LinearFVBCs]` and
  `[UserObjects]` blocks, routed by `isKokkosObject()`.
- Static polymorphism through `template <typename Derived> KOKKOS_FUNCTION` hooks. The
  registration macros detect which hooks are implemented, and loops that are not implemented are
  never launched.
- `Moose::Kokkos::Array` containers with shallow copy. Holders (`MeshHolder`, `SystemHolder`) and
  `ReferenceWrapper`/`Scalar` value binding.
- PETSc as the memory boundary: Kokkos vectors alias PETSc device arrays, and matrices are
  uploaded with COO. The COO preallocation is already done only once (`KokkosMatrix.K:35`).
- Linear FV BCs keep the precomputed affine **boundary relations**
  (`BoundaryRelation{coefficient, source}` for value and normal gradient). This is already the
  right GPU design: a BC is evaluated once per boundary face, and kernels read arrays.

### Changed relative to the current Kokkos linear FV path

| Today | Proposed | Why |
|---|---|---|
| Side geometry stored per element side, padded `(side, elem)`. Interior faces stored twice. No face ID. | Face SoA arrays plus `elem_faces(side, e)` | NS needs one mass flux per face. Halves face geometry memory. |
| Flux kernels loop over directed sides. Each side recomputes the face and atomically adds to its own row. | Face pass computes once; gather pass writes rows | Face work in NS is heavy (gradients, deviatoric terms, interpolation). Avoids FP64 atomics. Deterministic. |
| Matrix writes use `Matrix::find`, a binary search in the row. | Precomputed CSR value slots per face and per cell | Removes the search from the hot loop. |
| Boundary data indexed through `bc_face_index` from a per-kernel table | Boundary faces are contiguous per boundary ID, so `bc_face_index = f - range.begin` | Removes indirection. |
| No `gC`, `eCN`, skewness vector, ghost centroids, gradients, or solution reads in kernels | Precomputed face geometry; device gradients; cell value accessors | Required by every NS kernel. |
The existing `KokkosLinearFVDiffusion`, `KokkosLinearFVAdvection`, `KokkosLinearFVSource` and the
two functor BCs should be **migrated** to the new face API, not kept as a parallel path. Their
tests in `test/tests/kokkos/linearfvkernels` then become the regression baseline. This touches
recently merged upstream code (#32492), so it should be coordinated with its authors.

---

## 2. Mesh data structures (goals 3, 4, 5)

### 2.1 Index spaces

| Index | Range | Notes |
|---|---|---|
| `ContiguousElementID e` (the existing Kokkos ID) | `[0, n_local)` local cells grouped by block; `[n_local, n_local + n_ghost)` ghost layer | Unchanged from today's `Moose::Kokkos::Mesh`. |
| `FaceID f` (new) | `[0, n_faces)` | Every face with at least one local cell. Same set as CPU `_face_info`. |
| `LocalDofIndex` | Kokkos `System` local dof space | Unchanged. Linear FV variables have one dof per cell. |
| `ContiguousBoundaryID b` | existing | Used to group boundary faces. |

### 2.2 Element ("ElemInfo") SoA arrays

These are indexed by `e`, **including the ghost layer**. Ghost geometry is not cached today
(`KokkosMesh.K:304-330` sizes only `n_local`).

| Array | Type | Use |
|---|---|---|
| `elem_info[e]` | `ElementInfo{type, id, subdomain}` | existing |
| `elem_volume[e]` | `Real` (times coordinate factor) | existing, to be extended to ghosts |
| `elem_centroid[e]` | `Real3` | existing, to be extended to ghosts (needed for reconstructions and RZ) |
| `elem_faces(side, e)` | `Array2D<FaceID>` (LEFT layout, so `side` varies fastest) | **New. Goal 5.** Holds `invalid` for sides past the element type's count. |
| `elem_n_faces[e]` | `uint8_t` (or looked up from the element type) | loop bound |

**Face orientation inside the gather** is implicit: `sign = (face_elem[f] == e) ? +1 : -1`.
- That costs one extra load per face. The alternative is to pack the sign into the `FaceID` high
  bit, which saves the load and caps faces at 2^31 per rank with 32-bit indices.
- Recommendation: start with the implicit comparison and switch only if profiling shows it matters.

**Padded versus compressed sparse row (CSR) connectivity.**
- Padded `(max_sides, n_cells)` matches the existing `(side, elem)` layout, needs no offsets
  array, and gives a fixed stride for coalesced access.
- With first-order hex/tet/prism meshes, `max_sides <= 6`, so padding waste is small.
- A CSR `elem_face_offset/elem_face_ids` pair is only worth it for meshes that mix very different
  element types, such as polyhedra.
- Because elements are grouped by block, the padding can also be **per block** (one stride per
  block range) if mixed-topology blocks become common.

### 2.3 Face ("FaceInfo") SoA arrays

Indexed by `f` (goal 4):

| Array | Type | Notes |
|---|---|---|
| `face_elem[f]` | `ContiguousElementID` | Index into the element arrays (replaces `const ElemInfo *`). Always valid. |
| `face_neighbor[f]` | `ContiguousElementID` | Index or `invalid` on a domain boundary. May be a ghost ID. |
| `face_elem_side[f]`, `face_neighbor_side[f]` | `uint8_t` | Needed by BCs, the legacy `(side, elem)` access, and debugging. |
| `face_area[f]` | `Real` | Already multiplied by the face coordinate factor. |
| `face_normal[f]` | `Real3` | Outward from `elem`. |
| `face_gc[f]` | `Real` | Interpolation weight from `FaceInfo::gC()`. 1 on boundaries. |
| `face_dcn[f]` | `Real3` | `C_N - C_E`, or `x_f - C_E` on a boundary. |
| `face_inv_dn[f]` | `Real` | Precomputed `1 / abs(dCN . n)` or `1 / abs(dCN)`, depending on the non-orthogonal correction (NOC) setting. |
| `face_noc[f]` | `Real3`, optional | Precomputed `n - eCN / (n . eCN)`. Allocated only if a kernel uses NOC. |
| `face_centroid[f]` | `Real3`, optional | Only for position-dependent functions and wall distances. |
| `face_skew[f]` | `Real3`, optional | Only for skewness-corrected interpolation. |

`face_boundary[f]` is not stored; the boundary ID is implied by the face's range (2.4).

**Memory.**
- With every optional field allocated, a face costs about 140 bytes. Without them it costs about
  90 bytes.
- Precomputed fields trade memory for flops and divides. Each optional field should be allocated
  only on request, following the existing `setNeedsElementSideGeometry()` pattern.

**Building the arrays.** Build them on the host from the existing `MooseMesh` `FaceInfo`/`ElemInfo`.
Map each `FaceInfo`'s elem and neighbor to contiguous IDs through the extra-integer and ghost maps,
then copy to the device once.
- This reuses the CPU geometry exactly: `gC` from the face-plane intersection, vertex-average
  centroids, and coordinate factors. So the CPU and GPU see identical inputs, and CPU/GPU parity
  tests are meaningful.
- Rebuild in `Mesh::update()` whenever FV info changes.

### 2.4 Per-block organization (goal 3)

Faces are **sorted into contiguous groups**. Each group is described by a host and device range
table `FaceRange{begin, end, kind, block, other}`:

```
for each block s (in contiguous-subdomain order):
    [ interior faces of s, both cells local               ]  kind = INTERIOR
    [ interior faces of s, neighbor is a ghost cell        ]  kind = PROCESSOR
    [ boundary faces of s, grouped by boundary ID b       ]  kind = BOUNDARY, other = b
for each block pair (s1 < s2):
    [ faces between s1 and s2                             ]  kind = INTERFACE, other = s2
```

This gives a vectorized container per block **without nested arrays**.
- The "container" for block `s` is the slice of the flat SoA arrays given by its ranges.
- Element data is already contiguous per block.
- A block-restricted kernel's worklist is a short list of ranges. It needs no per-kernel
  `std::vector<FaceID>`, and the dispatch maps a thread ID to (range, offset) with a binary search
  over a handful of prefix sums.
- A boundary condition's faces are one range per boundary ID, so its relation arrays are indexed
  by `f - range.begin`.
- `faceType` (BOTH/ELEM/NEIGHBOR) becomes a property of the range for a given variable's block
  restriction. It is not a per-face nested vector. Interface faces are interior for a variable
  active on both blocks and boundary faces otherwise.
- Putting processor-boundary faces last allows **overlapping** the interior face pass with the
  ghost exchange later.
- *Within* a range, sort faces by `face_elem`, so face and gather passes walk cells roughly
  monotonically. That gives better cache and coalescing behavior.

**Alternative considered.** An `Array<BlockFVData>` with separate SoA allocations per block, in the
style of material property storage.
- Pro: per-block reallocation and per-block strides.
- Con: one dispatch per block or a nested indirection on device, more allocations, and harder
  cross-block (interface) faces.
- Recommendation: flat arrays plus ranges. See open question Q1.

### 2.5 Orientation and processor boundaries

- Use the CPU ownership rule (`Moose::FV::elemHasFaceInfo`): `elem` is the lower **global** ID or
  the finer side.
  - Ranks with different contiguous numberings then agree on the orientation of every shared face.
  - With identical ghosted inputs, the same arithmetic is performed on both ranks. So the mass
    flux on a processor-boundary face has the same value on both sides, and continuity holds
    without communicating face data. This mirrors what the CPU `RhieChowMassFlux` relies on today.
- Each face writes only the rows of its **local** cells (`elem < n_local`, `neighbor < n_local`).
  - Ghost rows are never touched. FV-only Kokkos systems can therefore drop ghost rows from the
    CSR and from the COO upload.
  - Today ghost rows are allocated, preallocated in COO, and uploaded as zeros.

---

## 3. Field data on device

| Field | Storage | Ghosted? | Producer |
|---|---|---|---|
| Solution, old states | Kokkos `Vector` aliasing PETSc (existing) | yes | solver |
| Cell gradients (per variable, `dim` components) | PETSc-backed system vectors exposed as Kokkos vector tags (same layout idea as the CPU `LinearFVGradientManager`) | yes, via PETSc ghost update | device Green-Gauss |
| `HbyA`, `Ainv` (cell, per component) | PETSc vectors, device type | yes | PETSc ops in `computeHbyA` |
| Cell properties (`rho`, `mu`, `cp`, `k`, `mu_eff`, ...) | `Array<Real>(n_local + n_ghost)` | computed on ghosts too | Kokkos property objects (section 6) |
| Face mass flux `mdot[f]` | `Array<Real>(n_faces)`, restartable | not needed (D5) | Kokkos RC |
| `HbyA_f[f]`, `rhoAinv_f[f]` | `Array<Real>` and `Array<Real3>` | not needed | Kokkos RC |
| Boundary relations | `Array<Real>` per BC (existing) | n/a | Kokkos BCs |

**Why gradients live in PETSc-backed vectors.**
- They need ghost values for interior faces next to processor boundaries.
- PETSc ghost updates on device vectors, with GPU-aware MPI, give that without a custom halo
  exchange, consistent with "PETSc is the memory boundary".
- This also matches the CPU storage, so it stays restart-compatible.
- **Risk:** libMesh ghosted vectors are created with `VecCreateGhost`. Device support for
  `-vec_type kokkos` with ghosted vectors must be verified (see Q4). The fallback is a Kokkos halo
  exchange that reuses the System's existing ghost communication lists.

**Face-field restart.**
- `mdot` is declared with `declareRestartableData<Moose::Kokkos::Array<Real>>` and wrapped in a
  `ReferenceWrapper`, the existing pattern in `test/src/kokkos/kernels/KokkosRestartDiffusion.K`.
- Like the CPU map keyed by a local face ID, this is partition-dependent.

---

## 4. Assembly design (goal 1)

### 4.1 Passes per linear system

```
1. BC pass         for each Kokkos BC: fill boundary relations          (per boundary range)
2. Face pass       for each flux kernel: accumulate face coefficients   (per face range)
3. Gather pass     per local cell: diag, rhs from faces + elemental     (per cell)
                   + equation relaxation, if requested
4. Upload          MatSetValuesCOO / vector aliasing                     (existing path)
```

**Face pass.** One thread per face. For an internal face `f` with cells `E`, `N`, a kernel returns
an `InternalFaceContribution`:

```cpp
struct InternalFaceContribution
{
  Real elem;          // a_EE: coefficient of phi_E in row E
  Real neighbor;      // a_EN: coefficient of phi_N in row E
  Real rhs_elem;      // b_E
  Real rhs_neighbor;  // b_N
};
```

This mirrors the CPU `computeElem/NeighborMatrixContribution` and
`computeElem/NeighborRightHandSideContribution`, with the CPU antisymmetric convention: row N is
`(-a_EE, -a_EN)`. CPU kernels can therefore be ported almost line by line. Kernels that are not
antisymmetric can also define `neighborRow` to give `(a_NE, a_NN)` explicitly.

The face pass accumulates into per-system face scratch, adding over all flux kernels of the
variable:

- `off_EN[f] += a_EN` is written straight into the CSR value slot `slot_EN[f]` (if `E` is local).
- `off_NE[f] += -a_EE`, or the explicit `a_NE`, is written into `slot_NE[f]` (if `N` is local).
- `diag_E[f]`, `diag_N[f]`, `rhs_E[f]`, `rhs_N[f]` go into scratch arrays of four `Real`s per face.

Every write targets a slot owned by exactly one face, and kernels run one after another. So **no
atomics are needed**.

**Gather pass.** One thread per local cell `e`:

```
diag = sum over sides of (sign > 0 ? diag_E[f] : diag_N[f]) + elemental kernels' matrix terms
rhs  = sum over sides of (sign > 0 ? rhs_E[f]  : rhs_N[f])  + elemental kernels' rhs terms
     + boundary contributions (face pass on boundary ranges writes diag_E/rhs_E the same way)
matrix(slot_diag[e]) = diag;  rhs[e] = rhs
```

Elemental kernels (sources, time derivative, Boussinesq, friction, pressure gradient) are
evaluated **inside** the gather through the same templated-hook mechanism. Either the system
gather calls each elemental kernel's hook, or each elemental kernel keeps its own cell pass with
non-atomic writes as today; the second is simpler for the dispatcher. Both are race-free because
each cell is one thread.

**Precomputed slots.** `slot_diag[e]`, `slot_EN[f]` and `slot_NE[f]` are CSR value indices
computed once from the sparsity pattern at setup. This replaces `Matrix::find` in the hot loop.

**Equation relaxation fused into the gather.**
- The gather already sees every off-diagonal of the row, so it can form
  `D* = max(|D|, sum |a_off|) / lambda` and add `(D* - D) * x_old` to the right-hand side directly.
- This removes the host `get_row` loop of `relaxMatrix` and two extra device passes.
- `RhieChowMassFlux::computeHbyA` must then read the **relaxed** diagonal. That is what the CPU
  does today, since it calls `MatGetDiagonal` after relaxation.

### 4.2 Alternatives considered

| Scheme | Face work | Atomics | Deterministic | Notes |
|---|---|---|---|---|
| Directed sides with atomics (current Kokkos path) | 2x | yes (diagonal, off-diagonal, rhs) | no | Simplest. Fine for cheap kernels such as constant-coefficient diffusion. |
| Face-centric scatter with atomics on diagonal and rhs | 1x | diagonal, rhs | no | One pass. FP64 atomics are fast on A100/H100/MI250 but slow elsewhere. |
| **Face pass + gather (proposed)** | 1x | none | yes | One extra pass, plus 4 `Real`s per face of scratch. |
| Cell-centric recompute (each cell evaluates its own faces) | 2x | none | yes | No scratch. Wins only when face work is trivial. |
| COO with duplicate entries summed by PETSc | 1x | none | yes | Very simple kernels. Pushes about 4·n_faces + n_cells values through PETSc's COO reduction on every assembly. Needs `VecSetValuesCOO` for the rhs. A candidate to benchmark against the proposed scheme. |

Bitwise reproducibility across runs is valuable when debugging SIMPLE convergence, which is a
second reason to avoid atomics.

### 4.3 Kernel API sketch

```cpp
class KokkosLinearWCNSFVMomentumFlux : public Moose::Kokkos::LinearFVFluxKernel
{
public:
  template <typename Derived>
  KOKKOS_FUNCTION InternalFaceContribution
  computeInternalFaceContribution(const FaceDatum & face) const;

  template <typename Derived>
  KOKKOS_FUNCTION BoundaryFaceContribution  // {diag, rhs} for the local side
  computeBoundaryFaceContribution(const FaceDatum & face, const BoundaryRelations & bc) const;
};
```

- `FaceDatum` is a value type built from `f` and the mesh. It gives the geometry accessors
  `area()`, `normal()`, `gC()`, `dCN()`, `invDn()`, `nocVector()`, the cells `elem()` and
  `neighbor()`, `isLocal(side)`, and cell accessors through the variable (`_u.elem(face)`,
  `_grad_u.neighbor(face)`).
- It replaces the CPU `setupFaceData` mutation. Per-face cached values become local variables
  inside the hook, so **there is no mutable per-face object state**.
- Matrix and right-hand side come from **one** hook, so face geometry, cell values and
  interpolation weights are read once.
  - The CPU `LinearFVComputationMode` (`RHS`, `Matrix`, `FullSystem`) stays available as a
    compile-time branch for RHS-only reassembly.
  - Today's Kokkos path has separate matrix and RHS loops, which read everything twice.
- **Shared face functions.** Physics math that two objects need is written once as a
  `KOKKOS_FUNCTION` free function, for example `pressureDiffusionFace(face, Ainv_f, grad_p)`. It is
  used by both the pressure diffusion kernel and `KokkosRhieChowMassFlux::computeFaceMassFlux`.
  Today RC instead calls a mutated kernel object's hooks.

### 4.4 Interpolation methods

- The `FVAdvectedInterpolationMethod` family is a set of pure functions of
  `(face geometry, phi_E, phi_N, grad_E, grad_N, flux)`.
- Port each as a `KOKKOS_FUNCTION` free function or a small trivially-copyable struct: upwind,
  average, harmonic, van Leer and minmod weight-based, MUSCL deferred correction.
- Kernels select the method **at compile time** through a template parameter or an `if constexpr`
  over an enum held by value. This avoids the device virtual calls Kokkos-MOOSE rules out.
- A runtime switch inside one kernel is acceptable for a small enum. GPU threads in one face range
  all take the same branch, so the branch is uniform and costs little.

### 4.5 Gradients on device

- Green-Gauss is a **cell gather** with no atomics:
  `grad[e] = (1/V_e) * sum over sides of sign * phi_f * n_f * A_f`.
- `phi_f` is the `gC`-weighted average on interior faces, or `coefficient * phi_C + source` from
  the variable's value relation on boundary ranges.
- The Venkatakrishnan limiter is a second cell pass.
- After computing, call the PETSc ghost update for the gradient vectors. Keep the CPU double
  buffering (compute all, then swap) so BCs see consistent previous gradients.

---

## 5. Rhie-Chow and the segregated solve on device

### 5.1 `KokkosRhieChowMassFlux`

It is a `GeneralUserObject`-style Kokkos object driven by the executioner (`execute_on = NONE`).

| CPU step | Device implementation |
|---|---|
| `computeHbyA` | PETSc device ops: `MatGetDiagonal`, `MatMult`, `VecPointwiseMult`, `VecAXPY`. The pressure-gradient removal is a cell kernel or `VecPointwiseMult`. SIMPLEC row sums come from the gather pass (4.1) or `MatGetRowSum`. |
| `populateCouplingFunctors` | Face pass over the flow ranges writes `HbyA_f[f]` and `rhoAinv_f[f]`. Boundary ranges use the velocity BC relations. |
| `computeFaceMassFlux` | Face pass: `mdot[f] = -HbyA_f[f] + pressureDiffusionFace(...)`, with the same shared function as the pressure kernel. |
| `updateCellVelocity` | Cell pass, or PETSc vector ops, using the device pressure gradient. |
| `getMassFlux(fi)` | `mdot[f]` read directly by kernels holding a `ReferenceWrapper` to the RC object (typed retrieval with `getKokkosUserObject<T>`). |
| `getVolumetricFaceFlux` | `mdot[f] / rho_f`, with `rho_f` interpolated from the cell property array. No string building or `isfinite` checks in the hot path. |

### 5.2 Executioner

Add a Kokkos-aware variant of `LinearAssemblySegregatedSolve`. It can be a mode of the existing
class rather than a new class (Q5).
- It keeps the CPU algorithm in section 5 of `ns_linearfv_design.md` unchanged and swaps the
  pieces listed above.
- All linear algebra stays in PETSc: KSP and preconditioner, normalization factor (two
  `MatMult`s), field relaxation (`VecAXPBY`), pressure pin.
- Users must select device matrix and vector types (`-vec_type kokkos`, `-mat_type aijkokkos`, or
  hypre for pressure). The executioner should **default** them for Kokkos systems, because a host
  type silently forces a full round trip per assembly.

### 5.3 Host-device traffic to eliminate

1. **Per-assembly synchronization.** `computeKokkosLinearSystem` syncs host to device and back on
   every assembly. With device PETSc vectors this becomes aliasing. Verify that FV-only systems
   take no host path in `Vector::close()`; D6 guarantees there are no ghost contributions.
2. **Per-assembly aux/user-object execution** in `FEProblemBase::computeLinearSystemTags`. Every
   object executed there must be a Kokkos object, or execution must be skipped for Kokkos-only
   systems.
3. **Host gradient computation** in `LinearSystem::computeGradients` for CPU variables. Kokkos
   variables must use the device path only.
4. **Solution ghost update after solves.** libMesh localizes `solution` to
   `current_local_solution`. This must be a device-aware PETSc scatter (Q4).
5. **Per-dispatch `Kokkos::fence()`.** For FV passes that are not followed by a host read, rely on
   stream ordering and fence once per linear system. That needs an opt-out in `Dispatcher`.

---

## 6. Properties and functors

CPU kernels evaluate virtual `Moose::Functor`s at `ElemArg`/`FaceArg`/`StateArg`. On device:

- **Cell property arrays.** A property (`rho`, `mu`, `cp`, `k`, `mu_eff`, `rho * cp`, ...) is
  materialized once per evaluation into `Array<Real>(n_local + n_ghost)` by a Kokkos FV property
  object, plus old states when time derivatives need them. Face values are interpolated inside
  kernels using `gC` (average) or harmonic averaging, matching the CPU `makeCDFace` path.
- **Sources supported, in order of effort:**
  1. constants (`Scalar`);
  2. `KokkosParsedFunction` of `(t, x)`;
  3. Kokkos parsed expressions of other cell fields (the existing RPN evaluator);
  4. ports of fluid-property user objects (a separate, larger effort).
- **Aux linear variables** (`T` from `h`, `mu_t`) need Kokkos aux kernels on linear FV variables:
  one dof per cell, so a cell pass. Check whether the existing Kokkos AuxKernel path handles
  `MooseLinearVariableFV`; it was written for FE variables.
- Properties evaluated on faces with a side (`functorFaceArg`, CHT) become one-sided reads of the
  cell array on the side given by the range kind.

---

## 7. Mapping the NS objects

| CPU object | Device notes | Phase |
|---|---|---|
| `LinearWCNSFVMomentumFlux` | Face pass. Reads `mdot[f]`, `mu` cells, own and all-component gradients (deviatoric terms), NOC vector. One hook for matrix and RHS. | 2 |
| `LinearFVPressureCorrectionDiffusion` | Face pass. `rhoAinv_f[f]` as the diffusion tensor; shares `pressureDiffusionFace` with RC. | 2 |
| `LinearFVDivergence` | Face pass, RHS only: `+-HbyA_f[f] * A`. | 2 |
| `LinearFVMomentumPressure` | Elemental: `-grad p * V` from the device gradient. | 2 |
| Framework BCs (Dirichlet, outflow, extrapolated) and NS BCs (extrapolated pressure, velocity and pressure symmetry, normal-velocity Dirichlet) | Boundary relations on boundary ranges. Velocity symmetry needs all components' cell values. Pressure symmetry reads `HbyA_f`. | 2 |
| `LinearFVTimeDerivative` | Elemental, with old property states. | 3 |
| `LinearFVEnergyAdvection`, `LinearFVScalarAdvection`, `LinearFVDiffusion` | Face pass with the shared interpolation functions. | 3 |
| Boussinesq, buoyancy, friction, volumetric heat transfer | Elemental. | 3 |
| `LinearFVPressureFluxBC` | Boundary relation reading `HbyA_f`, `rhoAinv_f`, `rho`. | 3 |
| Turbulence (TKE/TKED sources, turbulent advection and diffusion, wall-function BC, `mu_t` aux) | Wall data becomes cell and face arrays built at setup: a wall-bounded flag per cell plus a per-cell range of wall faces. Iterative y+ becomes a device Newton loop with a fixed maximum iteration count. The "skip wall rows" logic goes in the gather. | 4 |
| CHT (`CHTHandler`, Robin/convective BCs) | Interface face ranges with face arrays for flux and temperature. | 4 |
| Two-phase drift flux, slip | Verify the sign issue noted in `ns_linearfv_design.md` first. | 4 |
| `FVReconstructedPressureGradient` | Per-cell small least-squares solve with a fixed-size stack system: `dim x dim`, Cholesky via `Utils::choleskySolve`. | 4 |

---

## 8. Scalability (goal 1)

- **Inter-node.** There is no off-process assembly, and face data is not communicated.
  Communication per SIMPLE iteration is limited to:
  - ghost updates of the solution, gradients, `HbyA` and `Ainv`;
  - KSP communication;
  - global reductions for norms and convergence.

  This is the minimum for a cell-centered FV method. Use GPU-aware MPI through PETSc.
- **Intra-GPU.**
  - All hot loops are flat range loops over faces or cells with SoA, mostly coalesced access.
  - Face ranges are sorted by `face_elem`, which keeps cell gathers local.
  - There are no atomics and no searches.
  - Matrix and RHS are formed in one face pass per kernel.
- **Kernel launches.** For moderate meshes per GPU, time is bounded by launch latency (tens of
  kernels per equation per iteration). Mitigations, in order:
  1. drop per-dispatch fences (5.3);
  2. single combined matrix and RHS hooks (4.3);
  3. later, fuse all flux kernels of a variable into one face pass. This needs a composition
     mechanism in the dispatcher; Kokkos-MOOSE has none today.
- **Momentum components.** The u, v and w matrices share all face coefficients except
  boundary-condition rows and the deviatoric RHS.
  - A later optimization is to assemble the face pass once and copy values across the three
    systems (same nonzero pattern), patching boundary rows.
  - The systems stay separate MOOSE systems, so the Rhie-Chow diagonal per component is
    unaffected.
- **Solvers.**
  - Pressure: algebraic multigrid on device (hypre BoomerAMG or PETSc GAMG). Pressure is
    typically 50-80% of SIMPLE time, so this dominates scaling.
  - Momentum and scalars: GMRES or BiCGStab with Jacobi or block-Jacobi ILU.
  - Keep the preconditioner-reuse options.
- **Load balance.** Cells per rank decide the work, so the usual partitioners apply. Faces per
  rank follow the cells.
- **Precision and indices.** FP64 throughout. 32-bit indices are sufficient per rank.

---

## 9. Phased plan

| Phase | Content | Exit criterion |
|---|---|---|
| 0 | Face SoA, `elem_faces`, ghost cell geometry, face ranges, matrix slots, `FaceDatum`. Migrate the existing Kokkos linear FV kernels and BCs to the face pass + gather. | `test/tests/kokkos/linearfvkernels` and `fe-fv-interplay` give unchanged results. Atomics removed. |
| 1 | Device cell value accessors, Green-Gauss with limiter, NOC, `gC` interpolation, interpolation-method functions, cell property arrays (constant and parsed). | Kokkos `LinearFVDiffusion` with NOC and limited advection match the CPU to round-off on skewed meshes. |
| 2 | `KokkosRhieChowMassFlux`, momentum flux, pressure diffusion and divergence, momentum pressure, flow BCs, Kokkos-aware SIMPLE with fused relaxation and device PETSc defaults. | The `finite_volume/ins/{lid-driven,channel-flow,mms}/linear-segregated` inputs under `modules/navier_stokes/test/tests` match the CPU iteration history, and the MMS input keeps its convergence order. Profile with no host round trips. |
| 3 | PIMPLE, energy, scalars, Boussinesq and buoyancy, pressure-flux BC. | Corresponding CPU linear-segregated tests. |
| 4 | k-epsilon, CHT, two-phase, reconstructed pressure gradient, fluid-property ports. | Same. |

Use the existing `TIME_SECTION`s for each stage (copy, reinit, kernel, close) plus new ones for the
face pass, gather pass and gradients. Compare against the CPU at equal MPI rank counts and at one
rank per GPU.

---

## 10. Open questions

- **Q1. Per-block container shape.** Flat SoA arrays with per-block ranges (recommended), or a
  separately allocated SoA per block?
- **Q2. Face-array construction.** Build from `MooseMesh` `FaceInfo` (exact CPU parity,
  recommended), or natively in `Kokkos::Mesh` from libMesh? The native route avoids requiring
  `setupFiniteVolumeMeshData` but duplicates the geometry code.
- **Q3. Replace or coexist.** Migrate the existing upstream Kokkos linear FV `(side, elem)` path to
  the face design (recommended, needs upstream coordination), or add the face path beside it?
- **Q4. PETSc device ghosted vectors.** Does the PETSc/libMesh combination MOOSE pins support
  device ghost updates for `VecCreateGhost` vectors? This decides whether gradients and the
  solution ghost update stay in PETSc or need a Kokkos halo exchange.
- **Q5. Executioner shape.** A mode of `LinearAssemblySegregatedSolve`, or a separate Kokkos
  executioner class?
- **Q6. Target hardware.** NVIDIA only first, or AMD and Intel as well? This affects how
  attractive the atomic-based alternatives in 4.2 are as benchmarks.
