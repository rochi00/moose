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
  - "Per-face flux": `LinearFVAdvection` has a constant velocity only, so per-face flux storage
    comes with the Rhie-Chow mass flux in phase 2.
  - `bench.py parity --interp upwind average minmod vanleer muscl`: 240/240 match, matrices
    exactly. The steady bench inputs now start from a non-uniform initial condition and parity
    compares the first system (see the bench README for why the second does not work for Minmod
    and Van Leer).
  - Test: `kokkos/linearfvkernels/advection-interpolation`, one skewed-mesh advection-dominated
    problem per method (Van Leer and MUSCL with non-default parameters), against host golds.
