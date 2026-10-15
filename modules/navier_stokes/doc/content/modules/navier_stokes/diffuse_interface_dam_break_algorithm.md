# Diffuse-Interface Dam-Break Algorithm

This note documents the code path and numerical algorithm used by
`examples/dam_break/dam_break.i`. It focuses on the current reduced-pressure,
linear finite-volume, conservative diffuse-interface VOF path.

## Case Scope

The input `examples/dam_break/dam_break.i` configures a two-dimensional dam break:

- domain: `10a x 1.25a`
- mesh: `400 x 50`
- initial liquid column: `0 <= x < a`, `0 <= y < a`
- `a = 0.05715 m`
- liquid density: `rho_l = 998.19`
- gas density: `rho_g = 1.185`
- liquid dynamic viscosity: `mu_l = 1.0e-3`
- gas dynamic viscosity: `mu_g = 1.48e-5`
- gravity: `g = 9.81`
- timestep: `dt = 4.0e-4`
- VOF subcycles: `volume_fraction_subcycles = 1` (raised automatically if the alpha Courant limit demands)
- outer iterations: up to `num_iterations = 30`, stopping when the momentum and pressure residuals and the largest volume-fraction change are below `1e-6` (typically 7-10)
- relaxation: `momentum_equation_relaxation = 0.9`, `pressure_variable_relaxation = 0.7`
- alpha corrections: `n_alpha_corrections = 1`
- limiter iterations: `n_limiter_iterations = 6`

The solve uses four linear systems:

```text
u_system
v_system
pressure_system
alpha_system
```

The executioner is `ReducedPressurePIMPLE`. The key physics objects are:

- `WCNSLinearFVConservativeDiffuseInterfaceFlowPhysics`
- `WCNSLinearFVConservativeDiffuseInterfaceVOFPhysics`
- `ConservativeDiffuseInterfaceRhieChowMassFlux`
- `ConservativeDiffuseInterfaceVOFMULESCorrector`

## Governing Variables

The transported volume fraction is

\begin{equation}
\alpha \in [0, 1],
\end{equation}

where `alpha = 1` denotes liquid and `alpha = 0` denotes gas. Mixture properties
are generated as functors:

\begin{equation}
\rho(\alpha) = \alpha \rho_l + (1 - \alpha)\rho_g,
\end{equation}

and

\begin{equation}
\mu(\alpha) = \alpha \mu_l + (1 - \alpha)\mu_g.
\end{equation}

The flow solve uses reduced pressure, commonly written as

\begin{equation}
p_{rgh} = p - \rho \mathbf{g}\cdot(\mathbf{x} - \mathbf{x}_{ref}).
\end{equation}

In the code, the pressure variable is named `pressure`, but the diffuse-interface
flow physics requires the dynamic/reduced-pressure path.

## High-Level Algorithm

For each timestep, the reduced-pressure PIMPLE solve keeps this ordering:

```text
begin timestep
  initialize reduced-pressure / VOF transport state
  for each outer PIMPLE correction:
    advance outer-iteration histories
    solve alpha with VOF subcycling
    publish alpha-consistent rhoPhi
    solve momentum predictor
    solve pressure correction
    publish pressure-corrected velocity and face flux
  commit accepted transport history
end timestep
```

The important correctness property is that alpha is solved inside the outer
PIMPLE loop, before the momentum-pressure correction. This ensures the current
outer iteration sees updated `alpha`, `rho_mixture`, `mu_mixture`,
`alpha_phi_limited`, and `rho_phi`.

## Object Construction Path

The input selects these physics blocks:

```text
[Physics/NavierStokes/ConservativeDiffuseInterfaceFlowSegregated/flow]
[Physics/NavierStokes/ConservativeDiffuseInterfaceVOFSegregated/vof]
```

The flow physics enforces the conservative diffuse-interface Rhie-Chow object:

```text
WCNSLinearFVConservativeDiffuseInterfaceFlowPhysics
  -> rhieChowUserObjectType()
  -> "ConservativeDiffuseInterfaceRhieChowMassFlux"
```

The VOF physics creates the MULES corrector:

```text
WCNSLinearFVConservativeDiffuseInterfaceVOFPhysics::addUserObjects()
  -> addUserObject("ConservativeDiffuseInterfaceVOFMULESCorrector", "vof_mules", params)
```

The VOF corrector receives:

```text
system_name = alpha_system
variable = alpha
face_flux = vof_transport_phi
compression_factor = c_alpha
interface_normal = flow_interface_unit_normal_face
n_alpha_corrections = 1
n_limiter_iterations = 6
liquid_density = rho_l
gas_density = rho_g
```

The VOF physics also creates the matrix advection kernel for alpha:

```text
WCNSLinearFVConservativeDiffuseInterfaceVOFPhysics::addScalarAdvectionKernels()
  -> addLinearFVScalarAdvectionKernel(..., face_flux = "vof_transport_phi")
```

The top boundary gets an inlet-outlet scalar boundary condition, also driven by
`vof_transport_phi`.

## Published Face Fluxes

There are three central face-flux functors:

```text
corrected_face_phi
  pressure-corrected volumetric transport flux from Rhie-Chow

vof_transport_phi
  frozen volumetric flux used by the alpha equation during VOF transport

rho_phi
  alpha-consistent integrated mass flux published by the VOF corrector
```

The Rhie-Chow object owns `corrected_face_phi` and `vof_transport_phi`.
The VOF corrector owns `alpha_phi_limited` and `rho_phi`.

The VOF corrector computes the integrated alpha flux

\begin{equation}
\Phi_{\alpha,f}
\end{equation}

and publishes the integrated density flux

\begin{equation}
\Phi_{\rho,f}
  = \rho_g \Phi_f + (\rho_l - \rho_g)\Phi_{\alpha,f},
\end{equation}

where:

- `phi_f` is the volumetric face flux density from `vof_transport_phi`
- `A_f` is the face measure
- `Phi_f = phi_f A_f` is the integrated volumetric flux
- `Phi_{\alpha,f}` is the limited integrated alpha flux

The flow physics advertises `rho_phi` as the momentum mass-flux functor. The
diffuse-interface Rhie-Chow object therefore uses the VOF-owned `rho_phi` instead
of reconstructing mass flux as a simple face-density times volumetric-flux
product.

## Timestep Entry And Startup Projection

`ReducedPressurePIMPLE` installs `ReducedPressurePIMPLESolve` as the inner solve.
At the beginning of the segregated solve loop:

```text
LinearAssemblySegregatedSolve::solve()
  -> ReducedPressurePIMPLESolve::initializeSolveLoop()
```

The reduced-pressure initialization path is:

```text
ReducedPressurePIMPLESolve::initializeSolveLoop()
  -> resetVOFTransportStateForNewSolve()
     -> ConservativeDiffuseInterfaceRhieChowMassFlux::clearVOFTransportState()
  -> initializeConsistentStartupState()
     -> synchronizeSystemState(momentum systems)
     -> synchronizeSystemState(pressure system)
     -> synchronizeSystemState(alpha system)
     -> _problem.execute(EXEC_NONLINEAR)
  -> initializeStartupPressureField()
     -> assembleMomentumPredictorWithoutSolve()
     -> RhieChowMassFlux::initFaceMassFlux()
     -> correctStartupContinuityOnce() repeated startup_flux_corrections times
  -> commitAcceptedVOFTransportHistoryIfNeeded()
```

The startup projection is intentionally pressure-field preserving:

```text
ReducedPressurePIMPLESolve::correctStartupContinuityOnce()
  -> PressureStateGuard saves pressure vectors
  -> DiffuseInterfaceStartupProjectionGuard suppresses startup-only source terms
  -> preparePressureCorrectorState(true)
  -> applyPressureCorrectionStage(recompute_face_mass_flux = true,
                                  publish_pressure_corrected_state = false)
  -> restore pressure vectors
  -> pressure_system.computeGradients()
```

The projection repairs face fluxes before the first real outer correction, but
restores the user/equilibrium initial reduced-pressure field.

## Outer PIMPLE Loop

The base segregated driver owns the outer loop:

```text
LinearAssemblySegregatedSolve::solve()
  -> initializeSolveLoop()
  -> while outer_iteration < num_iterations:
       preMomentumPressureIteration()
       solveMomentumPredictor()
       correctVelocity()
```

For this dam-break case:

```text
num_iterations = 30
num_piso_iterations = 0
num_pressure_nonorthogonal_correctors = 0
momentum_absolute_tolerance = 1e-6
pressure_absolute_tolerance = 1e-6
active_scalar_absolute_tolerance = '1e-6'
```

So each timestep repeats the alpha update, the momentum predictor and the
pressure-correction stage until the normalized momentum and pressure residuals
and the largest change of the volume fraction since the previous outer
iteration are all below the tolerances, typically after 7-10 iterations; the
count is a cap. The volume fraction is judged by its change rather than by the
residual of its linear system, because the limited field never satisfies the
upwind system it was solved from, and the MULES limiter of a face may only
tighten from one outer iteration to the next, which keeps the loop from
cycling between two limiter states.

## Alpha Solve Ordering

The alpha solve is entered before the momentum predictor:

```text
ReducedPressurePIMPLESolve::preMomentumPressureIteration()
  -> advanceOuterIterationHistories()
  -> solveVolumeFractionBeforeFlowCorrection()
```

The VOF preparation and solve path is:

```text
ReducedPressurePIMPLESolve::solveVolumeFractionBeforeFlowCorrection()
  -> prepareVOFTransportStateForOuterIteration()
     -> clearVOFTransportState()
     -> freezeVOFTransportState(use_previous_timestep_transport_flux)
  -> _problem.execute(EXEC_NONLINEAR)
  -> solveVolumeFractionSystems()
  -> adoptPublishedVOFTransportState()
  -> _problem.execute(EXEC_NONLINEAR)
  -> storeActiveScalarResiduals()
```

`freezeVOFTransportState()` copies either the previous timestep corrected face
flux or the current corrected face flux into `vof_transport_phi`. On the first
timestep and first outer iteration, the previous timestep flux is used. Later
outer states use the current corrected flux.

## Alpha Subcycling

The alpha solve can subdivide the global timestep:

\begin{equation}
N_s = \max\left(
  N_{configured},
  \left\lceil \frac{Co_\alpha}{Co_{\alpha,max}} \right\rceil
\right),
\end{equation}

with

\begin{equation}
\Delta t_s = \frac{\Delta t}{N_s}.
\end{equation}

The executioner temporarily mutates problem time for each subcycle:

\begin{equation}
t_{old,s} = t^n + s\Delta t_s,
\end{equation}

\begin{equation}
t_s = t^n + (s + 1)\Delta t_s.
\end{equation}

The call path is:

```text
ReducedPressurePIMPLESolve::solveVolumeFractionSystems()
  -> computeVolumeFractionSubcycles()
  -> solveOneVolumeFractionSystem()
     -> system.saveOldSolutions()
     -> setPreviousNewtonToCurrent()
     -> corrector.resetSubcycleFluxes()
     -> for each subcycle:
          runOneVolumeFractionSubcycle()
     -> system.restoreOldSolutions()
     -> setPreviousNewtonToCurrent()
  -> restore global problem time
  -> finalizeVolumeFractionTransportState()
```

Each subcycle performs:

```text
ReducedPressurePIMPLESolve::runOneVolumeFractionSubcycle()
  -> setProblemSubcycleTime()
  -> advanceVolumeFractionSubcycleOldState() for subcycle > 0
  -> _problem.execute(EXEC_NONLINEAR)
  -> solveAdvectedSystem(alpha_system, ...)
  -> alpha_system.computeGradients()
  -> ConservativeDiffuseInterfaceVOFMULESCorrector::applyCorrection()
```

The matrix solve is donor/upwind alpha transport. The explicit correction is
applied afterward by the MULES corrector.

The true timestep-old alpha is restored after subcycling. This matters because
transient terms and stateful objects must still see the real old-time state, not
the intermediate subcycle state.

## Donor Alpha Equation

The base matrix alpha solve represents the conservative transport equation

\begin{equation}
\frac{V_P}{\Delta t_s}(\alpha_P^* - \alpha_P^n)
+ \sum_f \Phi_{f}^{upwind}\alpha_{upwind} = 0,
\end{equation}

where:

- `V_P` is the cell volume
- `alpha_P^*` is the donor/upwind solution after the linear solve
- `Phi_f` is the integrated volumetric face flux from `vof_transport_phi`
- the matrix interpolation is upwind

The solve is performed by:

```text
LinearAssemblySegregatedSolve::solveAdvectedSystem()
  -> computeLinearSystemSys()
  -> relaxMatrix()
  -> relaxRightHandSide()
  -> linear_solver.solve()
  -> system.setSolution()
```

In the reduced-pressure VOF path, the inherited scalar lower limiter is disabled
for this solve. Boundedness is enforced through the limited face-flux correction.

## MULES Correction Target Flux

The corrector starts from the donor flux:

\begin{equation}
\Phi_{\alpha,f}^{donor}
  = \Phi_f \alpha_{upwind}.
\end{equation}

For internal faces it computes a high-order face value using the Van Leer limited
face interpolation:

\begin{equation}
\alpha_f^{HO} = \alpha_f^{VanLeer}.
\end{equation}

For boundary faces, the high-order value is selected from the active boundary
condition. In the dam-break input the top patch uses an inlet-outlet alpha
condition; wall patches are treated as closed for alpha flux.

For internal faces with nonzero transport flux, the compression term modifies
the high-order value:

\begin{equation}
\alpha_f^{comp}
  = bound\left[
      \alpha_f^{HO}
      + c_\alpha \operatorname{sign}(\Phi_f)
        \alpha_f^{lin}(1 - \alpha_f^{lin})
        (\mathbf{n}_{\alpha,f}\cdot\mathbf{n}_f)
    \right],
\end{equation}

where

\begin{equation}
\alpha_f^{lin}
  = bound\left[g_C\alpha_P + (1 - g_C)\alpha_N\right].
\end{equation}

The target alpha flux is then

\begin{equation}
\Phi_{\alpha,f}^{target}
  = \Phi_f \alpha_f^{comp}.
\end{equation}

The raw explicit correction is

\begin{equation}
\Phi_{\alpha,f}^{corr}
  = \Phi_{\alpha,f}^{target}
    - \Phi_{\alpha,f}^{working}.
\end{equation}

On the first correction sweep, the working flux is the donor flux. Later sweeps
use the previously accepted working flux.

The corresponding call stack is:

```text
ConservativeDiffuseInterfaceVOFMULESCorrector::applyCorrection()
  -> collectFaceCorrectionData()
     -> buildFaceCorrectionData()
        -> faceTransportData()
        -> donorFlux()
        -> highOrderFaceValue()
           -> sharedVanLeerFaceValue() for internal faces
        -> compression calculation
```

## Limiter Math

For each cell, the corrector builds local admissible bounds from neighboring and
boundary alpha values, then limits the correction fluxes so the explicit update
does not leave those bounds.

The finite-volume explicit correction has the form

\begin{equation}
\alpha_P^{new}
  = \alpha_P^*
    - \frac{\Delta t_s}{V_P}
      \sum_{f \in \partial P} \lambda_f \Phi_{\alpha,f}^{corr},
\end{equation}

with face-orientation signs handled by the code through element/neighbor
updates. The limiter coefficient satisfies

\begin{equation}
0 \le \lambda_f \le 1.
\end{equation}

The available positive and negative correction budgets are built as

\begin{equation}
\psi_P^+
  = \frac{V_P}{\Delta t_s}\max(0, \alpha_{P,max} - \alpha_P^*),
\end{equation}

and

\begin{equation}
\psi_P^-
  = \frac{V_P}{\Delta t_s}\max(0, \alpha_P^* - \alpha_{P,min}).
\end{equation}

The code accumulates positive outflow and negative inflow correction sums for
each cell. It then iteratively tightens cell-side limiter estimates and applies
the minimum permitted value to each face. Only the process owning a cell sees
all of its faces, so the bounds, sums and limiters are computed for owned cells
only, and the limiters of cells next to a processor-partition face are sent to
the neighboring rank. Both adjacent ranks then accept the same face limiter,
and the result does not depend on the partitioning.

After limiting, the applied face correction is:

\begin{equation}
\Phi_{\alpha,f}^{limited}
  = w_k \lambda_f \Phi_{\alpha,f}^{corr},
\end{equation}

where

\begin{equation}
w_k =
\begin{cases}
1, & k = 0, \\
0.5, & k > 0.
\end{cases}
\end{equation}

For the dam-break input, `n_alpha_corrections = 1`, so only `w_0 = 1` is used.

The limiter/update call stack is:

```text
ConservativeDiffuseInterfaceVOFMULESCorrector::applyCorrection()
  -> build local bounds
  -> accumulate raw correction fluxes
  -> iterate n_limiter_iterations times:
       accumulate limited correction fluxes
       compute cell lambda_minus/lambda_plus of owned cells
       exchangePartitionCellLimiters()
       restrict accepted face lambda
  -> build limited_update vector
  -> current_local_solution.add(limited_update)
  -> publishFaceFluxes() on final correction
```

The published alpha flux is accumulated over subcycles with the same fraction
used for `rho_phi`:

\begin{equation}
\bar{\Phi}_{\alpha,f}
  = \sum_s \frac{\Delta t_s}{\Delta t}
      \Phi_{\alpha,f,s}^{limited}.
\end{equation}

## Alpha-Consistent Mass Flux

After correction, `publishFaceFluxes()` writes both:

```text
alpha_phi_limited
rho_phi
```

For each face and subcycle,

\begin{equation}
\Phi_{\rho,f,s}
  = \rho_g \Phi_{f,s}
    + (\rho_l - \rho_g)\Phi_{\alpha,f,s}^{limited}.
\end{equation}

The timestep-averaged published value is

\begin{equation}
\bar{\Phi}_{\rho,f}
  = \sum_s \frac{\Delta t_s}{\Delta t}\Phi_{\rho,f,s}.
\end{equation}

After all alpha systems finish, the executioner restores the global time and
calls:

```text
ReducedPressurePIMPLESolve::finalizeVolumeFractionTransportState()
  -> alpha_system.computeGradients()
  -> corrector.refreshPublishedRhoPhi()
```

`refreshPublishedRhoPhi()` recomputes `rho_phi` from the final published
`alpha_phi_limited` and current face state.

## Momentum Predictor

After alpha transport and nonlinear object execution, the base segregated solve
solves the momentum predictor:

```text
LinearAssemblySegregatedSolve::solve()
  -> solveMomentumPredictor()
```

The conservative diffuse-interface flow physics configures the momentum equation
to use `rho_phi` as its mass flux. The momentum predictor also receives an
explicit reduced-pressure diffuse-interface forcing through:

```text
ReducedPressurePIMPLESolve::addMomentumPredictorExplicitForcing()
  -> ConservativeDiffuseInterfaceRhieChowMassFlux::addMomentumPredictorExplicitForcing()
     -> reducedPressureMomentumPredictorForceDensity()
```

The force-density reconstruction includes pressure-gradient and, when not
suppressed, hydrostatic terms:

\begin{equation}
\mathbf{f}_{p,rgh}
  \approx -\nabla p_{rgh}
           - gh_f \nabla_n \rho.
\end{equation}

The startup projection guard suppresses these startup-only source pieces during
projection-only flux cleanup, but the normal solve enables the hydrostatic
predictor path.

## Pressure Correction

The PIMPLE pressure-correction path is:

```text
PIMPLESolve::correctVelocity()
  -> storePressurePreviousOuterIterationState()
  -> preparePressureCorrectorState(subtract_updated_pressure = true)
     -> RhieChowMassFlux::computeHbyA()
     -> ReducedPressurePIMPLESolve::postPreparePressureCorrectorState()
        -> updateAdditionalPressureFluxFunctors()
        -> updateVelocityBoundaryState()
        -> updatePressureBoundaryNormalGradients(false)
  -> applyPressureCorrectionStage(recompute_face_mass_flux = true,
                                  publish_pressure_corrected_state = true)
     -> solvePressureCorrector()
     -> postPressureCorrectorSolve()
        -> pressure_system.computeGradients()
        -> cachePressureEquationFlux()
     -> publishPressureCorrectedState()
        -> computeFaceMassFlux()
        -> relaxPressureFieldForNextPredictor()
        -> computeCellVelocity()
        -> postPublishPressureCorrectedState()
```

The diffuse-interface Rhie-Chow object constructs a predictor flux, pressure flux,
and pressure-corrected volumetric flux. Conceptually,

\begin{equation}
\phi_f^{corr}
  = \phi_f^{H/A} + \phi_f^p,
\end{equation}

where `phiHbyA` is the predictor contribution and `phi_p` is the pressure
equation contribution. The corrected volumetric flux is written to
`corrected_face_phi`.

The mass flux stored in the inherited Rhie-Chow face-mass field is computed from
the corrected volumetric flux and the interpolated face density:

\begin{equation}
\dot{m}_f
  = \rho_f \phi_f^{corr}.
\end{equation}

For downstream conservative diffuse-interface advection, however, the flow physics
uses the VOF-published `rho_phi` functor. That distinction is deliberate:
`rho_phi` is tied to the limited alpha flux and is not assumed to equal
`rho_f phi_f`.

## Transport History Commit

At timestep end:

```text
ReducedPressurePIMPLE::postTakeStep()
  -> ReducedPressurePIMPLESolve::commitAcceptedTimestepTransportHistory()
     -> ConservativeDiffuseInterfaceRhieChowMassFlux::
        commitAcceptedTimestepTransportHistory()
```

This copies the accepted `corrected_face_phi` into
`previous_timestep_corrected_face_phi`. The next timestep can then seed the
first alpha transport state from the accepted previous-timestep flux.

## Full Runtime Call Stack

The following stack omits framework scheduling details, but shows the relevant
dam-break path:

```text
ReducedPressurePIMPLE::ReducedPressurePIMPLE()
  -> _fixed_point_solve->setInnerSolve(ReducedPressurePIMPLESolve)

LinearAssemblySegregatedSolve::solve()
  -> ReducedPressurePIMPLESolve::initializeSolveLoop()
     -> resetVOFTransportStateForNewSolve()
     -> initializeConsistentStartupState()
     -> initializeStartupPressureField()
        -> assembleMomentumPredictorWithoutSolve()
        -> RhieChowMassFlux::initFaceMassFlux()
        -> correctStartupContinuityOnce()
           -> PIMPLESolve::preparePressureCorrectorState()
           -> PIMPLESolve::applyPressureCorrectionStage()
              -> solvePressureCorrector()
              -> postPressureCorrectorSolve()
              -> ConservativeDiffuseInterfaceRhieChowMassFlux::computeFaceMassFlux()
     -> commitAcceptedVOFTransportHistoryIfNeeded()

  -> outer iteration 1
     -> ReducedPressurePIMPLESolve::preMomentumPressureIteration()
        -> advanceOuterIterationHistories()
        -> solveVolumeFractionBeforeFlowCorrection()
           -> prepareVOFTransportStateForOuterIteration()
              -> ConservativeDiffuseInterfaceRhieChowMassFlux::clearVOFTransportState()
              -> ConservativeDiffuseInterfaceRhieChowMassFlux::freezeVOFTransportState()
           -> _problem.execute(EXEC_NONLINEAR)
           -> solveVolumeFractionSystems()
              -> solveOneVolumeFractionSystem()
                 -> ConservativeDiffuseInterfaceVOFMULESCorrector::resetSubcycleFluxes()
                 -> runOneVolumeFractionSubcycle()
                    -> LinearAssemblySegregatedSolve::solveAdvectedSystem()
                    -> ConservativeDiffuseInterfaceVOFMULESCorrector::applyCorrection()
                 -> runOneVolumeFractionSubcycle()
                    -> LinearAssemblySegregatedSolve::solveAdvectedSystem()
                    -> ConservativeDiffuseInterfaceVOFMULESCorrector::applyCorrection()
              -> finalizeVolumeFractionTransportState()
                 -> ConservativeDiffuseInterfaceVOFMULESCorrector::refreshPublishedRhoPhi()
           -> ConservativeDiffuseInterfaceRhieChowMassFlux::adoptPublishedVOFTransportState()
           -> _problem.execute(EXEC_NONLINEAR)

     -> solveMomentumPredictor()
        -> ReducedPressurePIMPLESolve::addMomentumPredictorExplicitForcing()
           -> ConservativeDiffuseInterfaceRhieChowMassFlux::
              addMomentumPredictorExplicitForcing()

     -> PIMPLESolve::correctVelocity()
        -> PIMPLESolve::preparePressureCorrectorState()
           -> RhieChowMassFlux::computeHbyA()
           -> ConservativeDiffuseInterfaceRhieChowMassFlux::
              updateAdditionalPressureFluxFunctors()
           -> RhieChowMassFlux::updateVelocityBoundaryState()
           -> RhieChowMassFlux::updatePressureBoundaryNormalGradients()
        -> PIMPLESolve::applyPressureCorrectionStage()
           -> solvePressureCorrector()
           -> PIMPLESolve::postPressureCorrectorSolve()
              -> ConservativeDiffuseInterfaceRhieChowMassFlux::cachePressureEquationFlux()
           -> PIMPLESolve::publishPressureCorrectedState()
              -> ConservativeDiffuseInterfaceRhieChowMassFlux::computeFaceMassFlux()
              -> PIMPLESolve::relaxPressureFieldForNextPredictor()
              -> ConservativeDiffuseInterfaceRhieChowMassFlux::computeCellVelocity()
              -> ReducedPressurePIMPLESolve::postPublishPressureCorrectedState()

ReducedPressurePIMPLE::postTakeStep()
  -> ReducedPressurePIMPLESolve::commitAcceptedTimestepTransportHistory()
```

With `volume_fraction_subcycles = 1` in the dam-break input, the
`runOneVolumeFractionSubcycle()` branch appears once per outer iteration;
the executioner raises the subcycle count when the alpha Courant number
would otherwise exceed its limit.

## Correctness Invariants

The implementation relies on these invariants:

1. Alpha transport runs before momentum-pressure coupling in every outer
   correction.
2. `vof_transport_phi` is frozen for alpha transport and is not changed in the
   middle of a subcycle.
3. The alpha matrix solve remains donor/upwind; bounded high-order behavior is
   supplied by the explicit limited correction.
4. `solutionOld()` for alpha is temporarily advanced between subcycles, then
   restored to the true timestep-old alpha.
5. `rho_phi` is published from the same limited alpha flux used to update alpha.
6. Startup continuity projection may repair face fluxes, but must restore the
   startup reduced-pressure field.
7. Pressure boundary normal gradients are refreshed after `HbyA` construction
   and before pressure equation assembly.
8. The accepted `corrected_face_phi` is committed at timestep end so the next
   timestep can seed alpha transport consistently. The boundary pressure normal
   gradients and face `1/A` are committed with it, so the first momentum
   predictor of the next timestep applies the same boundary pressure force the
   converged step ended with.

## Files In The Path

- `modules/navier_stokes/examples/dam_break/dam_break.i`
- `modules/navier_stokes/src/executioners/ReducedPressurePIMPLE.C`
- `modules/navier_stokes/src/executioners/ReducedPressurePIMPLESolve.C`
- `modules/navier_stokes/src/executioners/PIMPLESolve.C`
- `modules/navier_stokes/src/executioners/LinearAssemblySegregatedSolve.C`
- `modules/navier_stokes/src/physics/WCNSLinearFVConservativeDiffuseInterfaceFlowPhysics.C`
- `modules/navier_stokes/src/physics/WCNSLinearFVConservativeDiffuseInterfaceVOFPhysics.C`
- `modules/navier_stokes/src/userobjects/ConservativeDiffuseInterfaceVOFMULESCorrector.C`
- `modules/navier_stokes/src/userobjects/ConservativeDiffuseInterfaceRhieChowMassFlux.C`
- `modules/navier_stokes/src/userobjects/ConservativeDiffuseInterfaceRhieChowMassFlux.C`
