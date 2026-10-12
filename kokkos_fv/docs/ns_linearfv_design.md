# Navier-Stokes Linear FV: Existing Design

This describes the linear finite volume (linear FV) Navier-Stokes implementation in MOOSE `next`
at `410cc90c7e`. Paths are relative to the `moose/` repository root; `NS/` is
`modules/navier_stokes/` and `FW/` is `framework/`. Line numbers are approximate. Items marked
*(inferred)* are readings of the code rather than statements made in it.

The implementation has four layers:

1. **Framework mesh data**: `ElemInfo`/`FaceInfo` geometry and dof caches built by `MooseMesh`.
2. **Framework linear FV assembly**: `LinearSystem`, linear FV kernels/BCs, gradient manager,
   interpolation-method objects.
3. **NS objects**: momentum/pressure/energy/scalar/turbulence kernels and BCs, functor materials,
   `RhieChowMassFlux`, the Physics that create them.
4. **NS segregated executioners**: `SIMPLE`/`PIMPLE` driving `LinearAssemblySegregatedSolve`.

---

## 1. Mesh data: `ElemInfo` and `FaceInfo`

Built by `MooseMesh::setupFiniteVolumeMeshData()` (`FW/src/mesh/MooseMesh.C:4156`):
`buildFiniteVolumeInfo()` -> `computeFiniteVolumeCoords()` -> `cacheFaceInfoVariableOwnership()`
-> `cacheFVElementalDoFs()`. Rebuilt when the FV info is marked dirty (e.g. after AMR).

### Containers (`MooseMesh.h:1764-1782`)

| Member | Type | Contents |
|---|---|---|
| `_elem_to_elem_info` | `unordered_map<dof_id_type, ElemInfo>` | every active local **and ghost** element, keyed by libMesh element ID |
| `_elem_info` | `vector<const ElemInfo *>` | local elements only; order follows hash-map iteration |
| `_all_face_info` | `vector<FaceInfo>` | owning storage |
| `_face_info` | `vector<const FaceInfo *>` | faces with a local elem **or** neighbor |
| `_elem_side_to_face_info` | `unordered_map<pair<const Elem *, unsigned>, FaceInfo *>` | owning (elem, side) only |

### `ElemInfo` (`FW/include/mesh/ElemInfo.h`)

`const Elem *`, volume, centroid (vertex average), coordinate factor, and
`dofIndices()[sys][var]` (global PETSc dof index, `invalid_id` if the variable is inactive).

### `FaceInfo` (`FW/include/mesh/FaceInfo.h`)

- `const ElemInfo * _elem_info`, `const ElemInfo * _neighbor_info` (null on a mesh boundary).
- Dense local ID `_id` assigned sequentially while iterating active elements and sides
  (`MooseMesh.C:3876`). It is process-local.
- Geometry: area, face coordinate factor, normal (outward from *elem*), face centroid, `dCN`,
  `|dCN|`, `eCN`, interpolation weight `gC` (from the intersection of the C-N line with the face
  plane, so not 0.5 on skewed meshes), skewness-correction vector (computed on demand).
- Elem/neighbor side IDs, processor ID (that of *elem*), `std::set<BoundaryID>`.
- `faceType[sys][var]` in {BOTH, ELEM, NEIGHBOR, NEITHER}: on which sides the variable is active.
  A face where the variable lives on one side is treated as a boundary, which is how block
  boundaries (CHT, fluid/solid) are handled.

**Ownership rule** (`Moose::FV::elemHasFaceInfo`, `FW/src/utils/FVUtils.C:21-34`): the elem side
is the side with no neighbor, the finer side, or otherwise the lower element ID. A face is
*assembled* only by the rank owning its elem (`ownedFaceInfoBegin/End`), which then writes into
both rows; off-process rows are communicated by PETSc at `close()`.

**Implications for a GPU port:** pointer-linked AoS objects in hash-map order, a `std::set` and
nested `std::vector`s per face, and `unordered_map` lookups everywhere.

## 2. Framework linear FV assembly

### Systems and loops

`LinearSystem::computeLinearSystemInternal` (`FW/src/systems/LinearSystem.C:228-301`):
zero matrix/RHS -> optional `computeGradients()` -> Kokkos objects (if any) ->
`ComputeLinearFVElementalThread` over owned `ElemInfo` -> `ComputeLinearFVFaceThread` over owned
`FaceInfo` -> close. `FEProblemBase::computeLinearSystemTags` also executes multiapps, user
objects, and aux kernels **on every assembly**.

The CPU linear FV path is single-threaded in practice: `MooseLinearVariableFV` errors when
`numThreads() > 1` (`MooseLinearVariableFV.C:90`).

### Kernel API

- **`LinearFVElementalKernel`**: `computeMatrixContribution()`, `computeRightHandSideContribution()`;
  diagonal-only, inserted immediately with `matrix->add(dof, dof, v)`.
- **`LinearFVFluxKernel`**: `computeElemMatrixContribution()`,
  `computeNeighborMatrixContribution()`, `computeElemRightHandSideContribution()`,
  `computeNeighborRightHandSideContribution()`, `computeBoundaryMatrixContribution(bc)`,
  `computeBoundaryRHSContribution(bc)`, `computeBoundaryFlux(bc)`.
  - The face loop calls `setupFaceData(fi)` (kernels override it to cache per-face values),
    `setCurrentFaceArea(A * faceCoord)`, then add matrix/RHS.
  - **Internal face**: the flux is computed once and scattered antisymmetrically as a 2x2 block
    `[[a_E, a_N], [-a_E, -a_N]]` with `add_matrix`. Each row is written only if that side's
    subdomain is in the kernel's blocks. The RHS is `(rhs_E, rhs_N)`; kernels negate the neighbor
    RHS themselves.
  - **Boundary face** (variable on one side): BC found by
    `_var.getBoundaryCondition(*boundaryIDs().begin())` (map lookup; multiple boundary IDs on a
    face is an error), `bc->setupFaceData(fi, face_type)` (mutates the BC), diagonal-only add.
- **`LinearFVBoundaryCondition`**: `computeBoundaryValue()`, `computeBoundaryNormalGradient()`.
- **`LinearFVAdvectionDiffusionBC`**: affine relations
  - value: `phi_f = ValueMatrix * phi_C + ValueRHS` (used by advection),
  - gradient: flux `= k * (GradMatrix * phi_C - GradRHS)` (used by diffusion),
  - flags `includesMaterialPropertyMultiplier()`, `providesCompleteBoundaryFlux()`,
    `needsBoundaryNonorthogonalCorrection()`.
  Kernels `dynamic_cast` BCs to this type.

### Variables and functors

- Cell value: `getElemValue(elem_info, state)` reads the ghosted solution at the cached dof.
- Face value (functor evaluation): linear interpolation on internal faces; on boundary faces it
  calls `bc->setupFaceData(); bc->computeBoundaryValue()`, a hidden mutation of the BC.
- Kernels inherit `NonADFunctorInterface`; properties are `Moose::Functor<Real>` evaluated with
  `ElemArg`, `FaceArg` (`makeCDFace`, `singleSidedFaceArg`, `functorFaceArg`), and `StateArg`.

### Gradients (`LinearFVGradientManager`, `LinearFVGradientReader`)

- Variables call `requestCellGradients(method, oldest_state)` and get a reader.
- Storage per gradient *method*: `dim` ghosted `NumericVector`s (clones of the solution, so the
  gradient component of a variable sits at its own dof), with double buffering
  (`next_values` computed for all methods, then swapped) and named system vectors for old states.
- **Green-Gauss**: face pass over owned faces (`gC`-weighted linear interpolation; BC value or cell
  value on boundaries) scatter-adds `phi_f * n * A` into elem and neighbor via `add_vector`; volume
  pass divides by volume and applies the RZ correction. Optional Venkatakrishnan limiter pass.
  No skewness correction or iteration.
- Readers: `component(elem_info, d)` is a per-call `NumericVector` element read; `gradient(fi)`
  interpolates with `gC`.

### Interpolation-method objects (`FW/src/fvinterpolationmethods/`)

`FVAdvectedInterpolationMethod::advectedInterpolate(fi, phiE, phiN, gradE*, gradN*, flux)`
returns `{weights (wE, wN), rhs_face_value}`: `FVAdvectedUpwind`, `FVGeometricAverage`,
`FVHarmonicAverage`, `FVAdvectedVanLeerWeightBased`, `FVAdvectedMinmodWeightBased` (implicit
limited weights), `FVAdvectedMUSCLDeferredCorrection` (upwind matrix + deferred-correction RHS).
They are virtual MooseObjects but each is a pure function of face geometry, two cell values,
two gradients, and the flux.

## 3. NS kernels

### Momentum: `LinearWCNSFVMomentumFlux` (one per component)

- `setupFaceData`: mass flux `mdot_f = rc.getMassFlux(fi)` (hash-map lookup), advected
  interpolation weights, boundary normal factor.
- Internal matrix: `elem = (wE * mdot + mu_f / d) * A`, `neighbor = (wN * mdot - mu_f / d) * A`,
  with `mu_f` from the `mu` functor on a CD face; `d = |dCN . n|` with non-orthogonal correction,
  else `|dCN|`.
- Internal RHS: `+-(stress_rhs + rhs_face_value * mdot) * A`, where `stress_rhs` adds the
  non-orthogonal correction `mu_f (gC grad u_C + gN grad u_N) . (n - eCN / (n . eCN))` and, with
  `use_deviatoric_terms`, the explicit transpose/trace terms built from **all** velocity
  component gradients.
- Boundary: diffusion via the BC gradient relation (times `mu` unless the BC includes it),
  boundary non-orthogonal and deviatoric terms; advection via the BC value relation times `mdot`.
- Gradients requested: own (for NOC/deviatoric), advected-method gradient, all velocity
  components (deviatoric).

### Pressure

- `LinearFVPressureCorrectionDiffusion`: `LinearFVAnisotropicDiffusion` with the diffusion tensor
  being the `"Ainv"` face functor from `RhieChowMassFlux`. Matrix `(n o D) . n / d * A`; RHS
  holds tangential and non-orthogonal terms using pressure gradients.
- `LinearFVDivergence`: RHS-only `+-HbyA_f * A`, forced on boundaries.
- `LinearFVMomentumPressure` (elemental, RHS only): `-dp/dx_i * V` from the pressure gradient
  reader. With `FVReconstructedPressureGradient` it reads a relaxed coupling gradient at state 1.

### Other elemental kernels

Time derivative (`rho` factor, needs old functor states), Boussinesq, buoyancy, Darcy friction,
RZ viscous source, volumetric heat transfer, phase change, k-epsilon TKE/TKED sources. The TKE
kernels use per-element `std::map<const Elem *, ...>` wall data, allocate `std::vector`s per cell,
and run iterative y+ solves.

### Other flux kernels

`LinearFVEnergyAdvection` (legacy `InterpMethod` enum; constant `cp`), `LinearFVScalarAdvection`
(volumetric flux `mdot / rho_f`; optional AD slip-velocity functors),
`LinearFVTurbulentAdvection`/`LinearFVTurbulentDiffusion` (skip wall-bounded rows; write directly
to the system matrix, bypassing tags), framework `LinearFVDiffusion`, and
`LinearWCNSFV2PMomentumDriftFlux`.

### BCs

Framework Dirichlet/Neumann/Robin/extrapolated/outflow/symmetry, plus NS
`LinearFVExtrapolatedPressureBC`, `LinearFVNormalVelocityFunctorDirichletBC`,
`LinearFVVelocitySymmetryBC` (needs all velocity components), `LinearFVPressureSymmetryBC`
(cancels `HbyA`), `LinearFVPressureFluxBC` (prescribed mass flux, needs `HbyA`, `Ainv`, `rho`),
convective/Robin/Dirichlet CHT BCs, and `LinearFVTurbulentViscosityWallFunctionBC` (iterative
wall functions on the `mu_t` aux variable).

### Property sources

Constants, `NonADGeneralFunctorFluidProps` (fluid-property user objects, virtual calls),
`ADParsedFunctorMaterial`/`ParsedFunctorMaterial`, lambda-based functor materials
(`LinearFVEnthalpyFunctorMaterial`, mixture materials), aux linear variables (`T` from `h`,
`mu_t`), and the `RhieChowMassFlux` face maps.

## 4. `RhieChowMassFlux`

Executed only by the executioner (`execute_on = NONE`).

### Storage

| Data | Type |
|---|---|
| `_face_mass_flux` | `FaceCenteredMapFunctor<Real, unordered_map<dof_id_type, Real>>`, keyed by `FaceInfo::id()`, restartable |
| `_HbyA_flux` | same type, `rho_f * HbyA_f . n` |
| `_Ainv` | `FaceCenteredMapFunctor<RealVectorValue, ...>`, `rho * Ainv` per component |
| `_HbyA_raw[dim]`, `_Ainv_raw[dim]` | ghosted `NumericVector` clones of the momentum solution layout |
| `_grad_p_current[dim]` | snapshot of the pressure gradient used in the momentum predictor |
| `_cell_volumes` | `V * coordFactor` in pressure layout |
| `_flow_face_info` | `vector<const FaceInfo *>` of faces touching the flow blocks |

The face maps are keyed by a dense local face ID, so a plain face-indexed array would be a direct
replacement *(inferred: restart with different partitioning would also break today)*.

### Operations

- **`computeHbyA()`** (`RhieChowMassFlux.C:951-1151`), per component, on the assembled and
  relaxed momentum matrix `M`:
  - `A = diag(M)`, `H = M u - A u - b - grad p_snapshot * V` (removes the pressure source that
    `LinearFVMomentumPressure` put in `b`), `HbyA = H / A`, `Ainv = V / A`.
  - SIMPLEC (`pressure_projection_method = consistent`): row sums via host `get_row` loops.
- **`populateCouplingFunctors()`** (`:839-949`): serial host loop over flow faces. Internal faces
  interpolate `rho`, `HbyA`, `rho * Ainv` (average or harmonic); boundary faces use one-term
  values, with a Dirichlet-velocity special case that uses the boundary velocity.
- **`computeFaceMassFlux()`** (`:754-824`): serial host loop; for each face calls the pressure
  diffusion kernel's hooks (`setupFaceData` on the thread-0 kernel, face area 1) and the pressure
  BC's hooks to obtain the pressure flux, then `mdot_f = -HbyA_flux_f + p_flux`.
- **`updateCellVelocity()`**: `u_C = -HbyA_C - Ainv_C * grad p_C` (pure vector ops).
- **Queries**: `getMassFlux(fi)` (hash lookup); `getVolumetricFaceFlux(fi)` evaluates `rho` on a
  face, runs `isfinite` checks, and builds a `std::string` on every call.
- **Parallel**: no face-data communication. Every rank computes fluxes for all faces it sees using
  ghosted vectors *(inferred)*.

`FVReconstructedPressureGradient` (optional): per-cell least-squares projection of face fluxes to
a cell velocity (dense solve per cell), inverted to a relaxed coupling pressure gradient.

`CHTHandler`: per-interface-face loops evaluating BC/kernels and relaxing heat flux and interface
temperature stored in face maps.

## 5. Segregated executioners

All equations are separate `LinearSystem`s (one per velocity component, pressure, energy, solid
energy, each turbulence and scalar variable), each with its own PETSc matrix, RHS, and KSP. The
momentum KSP is shared across components and preconditioner reuse is configurable.

### One SIMPLE iteration (`LinearAssemblySegregatedSolve::solve`, `:743-987`)

1. First iteration: compute the pressure gradient. Snapshot it for momentum coupling.
2. **Momentum predictor** per component: assemble (with gradients) -> relax matrix/RHS ->
   normalization factor -> KSP solve.
3. **Pressure corrector**: `computeHbyA` -> `populateCouplingFunctors` -> assemble pressure ->
   optional pin -> KSP solve -> `computeFaceMassFlux` (unrelaxed pressure) -> field relaxation of
   pressure -> pressure gradient -> `updateCellVelocity`.
4. **Energy** (CHT fixed-point loop when enabled; fluid, radiation, solid).
5. Active scalars, then turbulence (with field relaxation and lower limiters).
6. `execute(EXEC_NONLINEAR)` for aux/user objects; convergence when every equation's normalized
   *initial* residual is below its tolerance.
7. Passive scalars in a separate loop after convergence.

**PIMPLE** reuses the same solve per time step and repeats step 3 `num_piso_iterations + 1`
times without reassembling momentum.

### Linear-algebra utilities (`SegregatedSolverUtils.C`)

- Equation relaxation: `D* = max(|D|, sum |offdiag|) / lambda`, `MatDiagonalSet`, and
  `b += (D* - D) x_old`. Row sums come from a host `get_row` loop.
- Field relaxation: `x = lambda x* + (1 - lambda) x_prev`.
- Normalization factor (OpenFOAM style): `|| |A x - A xbar| + |b - A xbar| ||_2`, two MatMults.
- Pressure pin: diagonal doubling plus RHS shift.

## 6. GPU readiness summary

| Already GPU-friendly | Needs redesign |
|---|---|
| Cell values and gradients (dof-indexed vectors) | `FaceInfo`/`ElemInfo` pointer graph and hash maps |
| Face geometry (normal, `dCN`, `gC`, area) | Face-indexed `unordered_map` storage in RC and CHT |
| Interpolation-method math (pure functions) | Virtual functors, fluid-property user objects, lambda materials |
| PETSc vector/matrix ops in RC, relaxation, normalization | Mutating `setupFaceData` used for evaluation |
| | Host `get_row` loops (relaxation, SIMPLEC) |
| | TKE wall maps, per-cell allocation, iterative wall functions |
| | Per-assembly aux/user-object execution |

### Points to verify before porting *(inferred)*

- `LinearWCNSFV2PMomentumDriftFlux`: the neighbor RHS does not appear to be negated relative to
  the elem RHS.
- `LinearWCNSFVMomentumFlux`: the boundary advection term has no outward-normal factor, unlike
  the scalar and energy kernels.
- `LinearFVTurbulentDiffusion`: the boundary non-orthogonal correction always uses the elem side.
- `force_boundary_execution` with no BC dereferences a null BC pointer (harmless only for
  `LinearFVDivergence`, which ignores it).
