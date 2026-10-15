# Phase 1 implementation plan: FV numerics and the rest of Group 0

**Started:** 2026-10-15. **Design:** `kokkos_ns_linearfv_detailed_design.md` sections 6, 7.2 and 14.

## Decisions (2026-10-15)

| Question | Decision |
|---|---|
| Scope | Phase 1 also takes the three Group 0 pieces the design placed later: reading other systems' variables, the time derivative, and the extra framework BCs (Robin, extrapolated, outflow, symmetry). Phase 2 can then start on SIMPLE directly. |
| BC parity | The Kokkos BCs match the host ones exactly: the Dirichlet distance is the normal projection `\|d_Cf . n\|`; the Neumann BC adds the tangential gradient term and takes a diffusive flux divided by the diffusion coefficient, as `LinearFVAdvectionDiffusionFunctorNeumannBC`. Kokkos gold files may change. |
| Gradients | G0C (design 6.2). |
| Kernel fusion, constant BC memoization | Later. |

## Steps

| Step | Content | Done when |
|---|---|---|
| 1. Reading variables | `LinearFVVariableInterface`: device handles to the cell values of a linear FV variable of any system, current, old or older; the assembly synchronizes the systems its objects read | A kernel reading a variable of another system reproduces the host |
| 2. Time derivative | `KokkosLinearFVTimeDerivative`, constant factor, implicit Euler and BDF2 | Matches `LinearFVTimeDerivative` for both integrators |
| 3. Gradients | G0C side table in `FVMesh`; Green-Gauss cell gather into the host gradient storage of `LinearFVGradientManager`; ghost update | Matches the host gradients; exact for linear fields on orthogonal meshes |
| 4. BC parity | Dirichlet distance; Neumann tangential term and flux semantics | `bench.py parity` boundary rows match on skewed, tet and prism meshes |
| 5. Non-orthogonal correction | Diffusion with NOC | Matches `LinearFVDiffusion` with `use_nonorthogonal_correction` on skewed meshes |
| 6. Advection | Per-face flux, interpolation methods, limiter | Matches `LinearFVAdvection` for each method |
| 7. Property arrays | Constant and parsed properties per cell, interpolated to faces | Used by diffusion and advection |
| 8. Framework BCs | Robin, extrapolated, outflow, symmetry | Each matches its host BC |

## Progress

- **Step 1 done.** `Moose::Kokkos::LinearFVVariableInterface` (`KokkosLinearFVVariableInterface.h/.K`)
  is a base of the Kokkos linear FV kernels and BCs. `getLinearFVCellValues(name, state)` returns a
  `LinearFVCellValues` handle (system, variable, solution tag), read on device with
  `handle(kokkosSystems(), elem)`, including one-layer ghost cells; it requests the solution state
  from the variable's system. `computeKokkosLinearSystem` synchronizes the union of the systems the
  objects read instead of only its own. Test: `kokkos/linearfvkernels/coupled-systems`, two linear
  systems, the second with a source read from the first, against the host objects.
- **Step 2 done.** `KokkosLinearFVTimeDerivative` with a constant `factor`, implicit Euler and BDF2
  (with its implicit Euler first step), using the BDF2 weights through a new const accessor
  `BDF2::weights()`. It replaces the test-only `KokkosLinearFVImplicitEulerTest`. Test:
  `kokkos/linearfvkernels/time-derivative`, both integrators, against host gold files.
- **Step 3 done (limiter later, with step 6).**
  - `FVMesh::requestSideData()` builds the G0C side table on request: per (owned cell, side),
    cell-fastest, the neighbor, the interpolation weight of the cell's own value and the outward
    area-weighted normal.
  - `KokkosFVGreenGaussGradient` is an `FVGradientMethod`, so the host `LinearFVGradientManager`
    keeps the storage, the old states and the ghost update, and host and Kokkos objects read the
    same gradients. Per variable it builds, once, the face value source of each (cell, side):
    internal, a Kokkos BC value relation, or the cell value where no BC applies (as
    `FVGreenGaussGradient`). Each owned cell sums its sides in a fixed order and writes its
    gradient into the owned entries of the component vectors; the RZ term is included.
  - Test: `kokkos/fvgradientmethods/green-gauss`, against an all-host run (the host Green-Gauss
    method does not see Kokkos BCs, so the reference must use host BCs too), serial, `-p 3` and
    distributed mesh. The mesh is non-uniform but orthogonal for now: on skewed meshes the Kokkos
    and host solutions still differ through the Dirichlet distance and the non-orthogonal
    correction (steps 4 and 5). Switch the test to a skewed mesh after them.
  - `LinearFVGradientComponentTestAux` (test app) outputs one gradient component, as the vector
    aux variables of `FunctorElementalGradientAux` are not supported with Kokkos FE.
- **Step 4 done.** The Kokkos and host linear systems now match to round-off in all 24 `bench.py
  parity` cases (hex, tet, prism; orthogonal and skewed; serial and 3 ranks), boundary rows
  included.
  - `KokkosLinearFVFunctorDirichletBC` uses `|d_Cf . n|` (`FaceDatum::cellToFaceDistance()`).
  - `KokkosLinearFVFunctorNeumannBC` takes the diffusive flux (`functor`) and a required
    `diffusion_coeff`, as the host BC; its boundary value adds the tangential gradient term, read
    through a new `LinearFVCellGradient` handle. BCs declare
    `includesMaterialPropertyMultiplier()`, which the diffusion kernel honours.
  - Gradient handles: `requestLinearFVCellGradient()` registers a member handle;
    `acquireGradients()`/`releaseGradients()` point it at the current gradient vectors (ghosted
    local form) around the dispatches of the assembly and of `KokkosFVGreenGaussGradient`.
  - Inputs with the Kokkos Neumann BC now give a flux and `diffusion_coeff`; the gold files did not
    change (orthogonal meshes, where `|d_Cf . n| = |d_Cf|` and the tangential term vanishes).
  - **Default gradient method (decided 2026-10-15): Kokkos automatically.** The host `green-gauss`
    method did not see Kokkos BCs and silently used cell values on their faces. Now
    `FVGreenGaussGradient` computes the gradients of the variables whose linear FV objects are all
    Kokkos objects with Kokkos (`LinearSystem::hasOnlyKokkosLinearFVObjects()`, set at the initial
    setup), and the others on the host as before; host-only problems are unaffected. The Kokkos
    computation is `Moose::Kokkos::GreenGaussGradient`, shared with `KokkosFVGreenGaussGradient`.
    Test: `kokkos/fvgradientmethods/host_method_name`; the bench inputs no longer select a method.
- **Step 5 done.** `KokkosLinearFVDiffusion` has `use_nonorthogonal_correction`, on by default as
  in `LinearFVDiffusion`, so Kokkos inputs on orthogonal meshes now compute gradients too.
  - Internal faces: the matrix coefficient divides by `|d_CN . n|`; the right-hand side gets
    `k (g_C grad_C + (1 - g_C) grad_N) . (n - d_CN / (n . d_CN)) A` on the elem row and its
    negative on the neighbor row. `FaceDatum` gained `faceGC()` and `neighborCentroid()`.
  - Boundary faces: BCs declare `needsBoundaryNonorthogonalCorrection()` (the Dirichlet BC does, as
    on the host); the kernel adds `k grad_C . (n - d_Cf / (n . d_Cf)) A` there.
  - Spherical coordinates: Green-Gauss gradients are not supported there (the host asserts in
    debug builds only; the Kokkos computation errors), so the `rspherical` MMS test disables the
    correction, which vanishes on its 1D mesh.
  - `bench.py parity` gained `--noc` (off and on by default) and compares the system of the second
    solve, where the gradients are nonzero in every cell: 48/48 match, matrices exactly, right-hand
    sides within 7e-16.
  - Test: `kokkos/fvgradientmethods/green-gauss` now runs on a skewed mesh with Dirichlet and
    Neumann BCs and the correction on, against an all-host gold; without the correction it fails.
- **Step 6 done.**
  - `KokkosLinearFVAdvection` takes the required `advected_interp_method_name`, as the host kernel,
    and supports the five framework methods: `FVAdvectedUpwind`, `FVGeometricAverage`,
    `FVAdvectedMinmodWeightBased`, `FVAdvectedVanLeerWeightBased` and
    `FVAdvectedMUSCLDeferredCorrection`. The constructor identifies the host method object and copies
    its parameters (`limit_to_linear`, `blending_factor`, `deferred_correction_factor`); the device
    code branches on the method, uniformly per kernel. Other methods are rejected with a parameter
    error. Existing Kokkos inputs gained an `FVAdvectedUpwind` block.
  - Methods that need gradients request them from the method's gradient method
    (`requestLinearFVCellGradient()` gained a method name), such as `green-gauss-venkatakrishnan`,
    the default of MUSCL.
  - Limiter on the device (decided 2026-10-15, option A): `FVGradientMethod` has a new protected
    virtual `applyLimiter()`, whose default is the former host limiting code of `computeGradient()`,
    so host-only problems and other gradient methods (`FVReconstructedPressureGradient`) run the
    same code. `FVGreenGaussGradient` (behind `MOOSE_KOKKOS_ENABLED`) overrides it to limit the
    Kokkos variables with `Moose::Kokkos::GreenGaussGradient::limit()` and the others with the
    base; `KokkosFVGreenGaussGradient` limits all on the device. The device limiter mirrors
    `ComputeLinearFVLimitedGradientThread` operation for operation over the G0C side table, with
    `Elem::hmin()` per cell precomputed. A separate kernel after the gradient kernel; fusing the two
    is possible since the limiter reads only the cell's own raw gradient (kernel fusion, later).
  - Tests: `kokkos/fvgradientmethods/venkatakrishnan` (limited gradients against a host gold; they
    differ from the unlimited ones); the MUSCL case of `advection-interpolation` now limits on the
    device and fails if the limiter is dropped.
- **Step 7 done.**
  - Survey (2026-10-15): Kokkos MOOSE has no reusable per-cell property storage for FV. Kokkos FE
    materials store values per quadrature point (or per element/subdomain with `constant_on`),
    owned elements only, computed inside the FE assembly loop; Kokkos functions are usable on GPU
    only through their concrete type. So linear FV gets its own property arrays.
  - `Moose::Kokkos::LinearFVPropertyInterface` (base of the Kokkos linear FV kernels and BCs):
    `requestLinearFVProperty(param, handle, use)` takes a `MooseFunctorName` that is a number, a
    `KokkosConstantFunction` or a `KokkosParsedFunction`. Constants live in the handle; functions
    fill device arrays over all owned and ghost cells and over the object's own faces, only when
    the time changes, before the object's dispatches. The device handle `LinearFVProperty` reads
    `cell(elem)` and `face(face)`.
  - Host evaluation points mirrored: cells at `ElemInfo::centroid()` (the vertex average, which is
    also `ElemArg::getPoint()`); faces at the face centroid, moved 1e-8 along the normal into the
    variable's cell on block-interface faces when the host uses `singleSidedFaceArg()` (diffusion
    coefficient, Neumann flux), not when it uses `functorFaceArg()` of a function (Dirichlet value).
    `coeff_interp_method` (`FVGeometricAverage`, `FVHarmonicAverage`) interpolates the cell values
    to the internal faces with the host formulas.
  - Converted: `KokkosLinearFVDiffusion` (`diffusion_coeff`, new `coeff_interp_method`),
    `KokkosLinearFVSource` (`source_density`, `scaling_factor`, default 1 as on the host),
    `KokkosLinearFVFunctorDirichletBC` (`functor`), `KokkosLinearFVFunctorNeumannBC` (`functor`,
    `diffusion_coeff` now default 1). The BC values are thus computed once per time, a first part
    of the static BC memoization. `KokkosLinearFVTimeDerivative` keeps a constant `factor`: a
    functor factor needs old states of the property (phase 3, with the density).
  - Found and fixed on the way: the boundary condition of a variable at a block interface is now
    selected from the boundary IDs of the variable's side, or else of the other side
    (`FVMesh::getFaceRangeVariableBoundaryIDs()`). The host takes the union of both sides and
    requires at most one ID, so a single interface sideset applies to the variables on both sides;
    Kokkos previously ignored it on the side without the sideset.
  - Test: `kokkos/linearfvkernels/properties`: numbers, a constant and parsed functions, harmonic
    and average interpolation, two block-restricted variables with a shared interface sideset and
    discontinuous interface functions that distinguish the sided and unsided evaluation, on a
    skewed mesh, against a host gold.
  - **Variable-based properties (decided 2026-10-15: additive parser change, exact host semantics,
    syntax mirroring `ParsedFunctorMaterial`).**
    - Shared Kokkos parser (`KokkosFunctionParser.h`, also used by Kokkos FE): a new `ARGUMENT`
      variable kind (`RPNBuilder::associateArgument()`) whose values the caller passes to
      `RPNEvaluator::eval()` through a new trailing `arguments` pointer, null by default. Existing
      call sites and results are unchanged (Kokkos FE function, material and parsed tests pass).
    - `KokkosParsedFunctorMaterial` (`[FunctorMaterials]`): the parameters of
      `ParsedFunctorMaterial` (`expression`, `functor_names`, `functor_symbols`, `property_name`),
      functors that are linear FV variables (solver or auxiliary) or numbers, `x y z t`, and `pi`,
      `e` with the six decimals the host passes. It defines no host functor.
    - Property parameters also accept a linear FV variable or a `KokkosParsedFunctorMaterial`
      `property_name`, resolved in `initialSetup()` (functor materials may be added after the
      kernels). Such properties are refilled before every assembly.
    - Host semantics mirrored (found by bisecting host/Kokkos differences): a functor material is
      evaluated on cells with the variables' cell values; on unsided internal and block-interface
      faces the host (`PiecewiseByBlockLambdaFunctor`) interpolates its cell values; on sided faces
      it evaluates the expression with the variables' face values, which are their linear
      interpolation, their boundary condition's boundary value, or their cell value
      (`MooseLinearVariableFV::evaluate()`). The boundary values come from the variable's Kokkos
      boundary conditions, recomputed before the fill with nested gradient acquisition
      (`acquireGradients()` now nests), followed by one fence. A property that needs itself through
      a boundary condition is an error. Variables with host boundary conditions are rejected.
    - The Green-Gauss computation now synchronizes the systems its boundary conditions read.
    - Ghosting: auxiliary linear FV variables had no ghost entries (empty send list), so their
      values on the ghost cells next to processor boundaries were unavailable to Kokkos. The Kokkos
      linear FV objects now declare a one-layer algebraic `ElementSideNeighborLayers` relationship
      manager, which MOOSE applies to every system. Note: the host reads such ghost values through
      a ghosted-vector lookup that is only checked in debug builds; its parallel results matched
      the serial ones here, but this deserves a look in a debug build.
    - Tests: `kokkos/linearfvkernels/properties/variables` (a functor material of a variable of
      another system with Dirichlet and Neumann BCs, an auxiliary variable, a number, `x` and `t`;
      and a variable directly) and `.../functor_material` (a functor material evaluated without a
      side on a block interface), against host golds, serial, `-p 3` and distributed mesh.
    - Not supported yet: Kokkos functions among `functor_names` (use `x y z t` in the expression),
      functor materials of functor materials, several functor materials defining one property on
      different blocks, old states of properties (needed by a functor time derivative factor).
  - "Per-face flux": `LinearFVAdvection` has a constant velocity only, so per-face flux storage
    comes with the Rhie-Chow mass flux in phase 2.
  - `bench.py parity --interp upwind average minmod vanleer muscl`: 240/240 match, matrices
    exactly. The steady bench inputs now start from a non-uniform initial condition and parity
    compares the first system (see the bench README for why the second does not work for Minmod
    and Van Leer).
  - Test: `kokkos/linearfvkernels/advection-interpolation`, one skewed-mesh advection-dominated
    problem per method (Van Leer and MUSCL with non-default parameters), against host golds.
