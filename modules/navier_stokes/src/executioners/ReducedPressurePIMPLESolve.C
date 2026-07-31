//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ReducedPressurePIMPLESolve.h"

#include "LinearSystem.h"
#include "ConservativeSharpInterfaceRhieChowMassFlux.h"
#include "ConservativeSharpInterfaceVOFMULESCorrector.h"
#include "MassConservedLevelSetCorrector.h"
#include "FEProblemBase.h"
#include "LinearFVLowMachDivergence.h"
#include "LinearFVLowMachEnthalpyAdvection.h"
#include "LinearFVLowMachMassAdvection.h"
#include "LinearFVLowMachSolidDrag.h"
#include "LinearFVTimeDerivative.h"
#include "LinearWCNSFVMomentumFlux.h"
#include "LowMachDivergenceSource.h"
#include "LowMachImplicitMassFlux.h"
#include "LowMachEnthalpyNewtonState.h"
#include "TheWarehouse.h"
#include "TimeIntegrator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
using namespace libMesh;

namespace
{
class ProblemTimeGuard
{
public:
  explicit ProblemTimeGuard(FEProblemBase & problem)
    : _problem(problem), _dt(problem.dt()), _time(problem.time()), _time_old(problem.timeOld())
  {
  }

  ~ProblemTimeGuard() { restore(); }

  void restore()
  {
    if (_restored)
      return;

    _problem.dt() = _dt;
    _problem.time() = _time;
    _problem.timeOld() = _time_old;
    _restored = true;
  }

private:
  FEProblemBase & _problem;
  const Real _dt;
  const Real _time;
  const Real _time_old;
  bool _restored = false;
};

class SharpInterfaceStartupProjectionGuard
{
public:
  explicit SharpInterfaceStartupProjectionGuard(ConservativeSharpInterfaceRhieChowMassFlux * rc)
    : _rc(rc),
      _saved_explicit_hydrostatic(rc ? rc->suppressExplicitHydrostaticPressureFlux() : false),
      _saved_startup_sources(rc ? rc->suppressStartupPressurePredictorFluxSources() : false)
  {
    if (_rc)
    {
      // Projection-only startup cleanup should solve only the continuity flux repair.
      _rc->setSuppressStartupPressurePredictorFluxSources(true);
      _rc->setSuppressExplicitHydrostaticPressureFlux(true);
    }
  }

  ~SharpInterfaceStartupProjectionGuard()
  {
    if (_rc)
    {
      _rc->setSuppressStartupPressurePredictorFluxSources(_saved_startup_sources);
      _rc->setSuppressExplicitHydrostaticPressureFlux(_saved_explicit_hydrostatic);
    }
  }

private:
  ConservativeSharpInterfaceRhieChowMassFlux * const _rc;
  const bool _saved_explicit_hydrostatic;
  const bool _saved_startup_sources;
};

class PressureStateGuard
{
public:
  explicit PressureStateGuard(LinearSystem & pressure_system)
    : _pressure_system(pressure_system),
      _linear_system(libMesh::cast_ref<LinearImplicitSystem &>(pressure_system.system())),
      _current_solution(*pressure_system.system().current_local_solution),
      _linear_solution(*_linear_system.solution),
      _saved_current(_current_solution.zero_clone()),
      _saved_linear(_linear_solution.zero_clone())
  {
    *_saved_current = _current_solution;
    _saved_current->close();

    *_saved_linear = _linear_solution;
    _saved_linear->close();

    if (auto * previous_newton = _pressure_system.solutionPreviousNewton())
    {
      _saved_previous_newton = previous_newton->zero_clone();
      *_saved_previous_newton = *previous_newton;
      _saved_previous_newton->close();
    }
  }

  ~PressureStateGuard()
  {
    if (!_restored)
      restore();
  }

  void restore()
  {
    _current_solution = *_saved_current;
    _current_solution.close();

    _linear_solution = *_saved_linear;
    _linear_solution.close();

    _pressure_system.setSolution(_current_solution);

    if (auto * previous_newton = _pressure_system.solutionPreviousNewton();
        previous_newton && _saved_previous_newton)
    {
      *previous_newton = *_saved_previous_newton;
      previous_newton->close();
    }

    _restored = true;
  }

private:
  LinearSystem & _pressure_system;
  LinearImplicitSystem & _linear_system;
  NumericVector<Number> & _current_solution;
  NumericVector<Number> & _linear_solution;
  std::unique_ptr<NumericVector<Number>> _saved_current;
  std::unique_ptr<NumericVector<Number>> _saved_linear;
  std::unique_ptr<NumericVector<Number>> _saved_previous_newton;
  bool _restored = false;
};
}

// Correctness invariants for this reduced-pressure sharp-interface solve:
//
// 1. The volume-fraction equation is solved inside each outer PIMPLE correction, before
//    momentum-pressure coupling, so rho, mu, rhoPhi, and VOF transport fluxes correspond to the
//    current outer state.
// 2. During VOF subcycling, solutionOld() is temporarily advanced between subcycles, but after
//    subcycling it must be restored to the true timestep-old alpha.
// 3. Startup continuity projection may repair face fluxes, but it must not overwrite the
//    user-supplied initial reduced-pressure field.
// 4. Problem time must be restored before publishing/finalizing VOF transport state.
// 5. Pressure boundary normal gradients must be refreshed after HbyA is computed and before the
//    pressure equation is assembled.
// 6. Explicit hydrostatic/startup pressure predictor flux suppressions are temporary startup
//    projection state and must always be restored.
// 7. Low-Mach enthalpy is converged and thermodynamically synchronized before momentum-pressure.
// 8. The temporary transported density is synchronized to the EOS only after the coupled solve is
//    accepted; the already-published common mass flux is not reconstructed during synchronization.

InputParameters
ReducedPressurePIMPLESolve::validParams()
{
  InputParameters params = PIMPLESolve::validParams();
  params.set<unsigned int>("num_iterations") = 1;
  params.setDocString(
      "num_iterations",
      "The number of outer PIMPLE corrections. For the transient reduced-pressure sharp-interface "
      "path this should remain a small outer-correction count, not a large SIMPLE-style "
      "momentum-pressure convergence loop. The reduced-pressure executioner performs this many "
      "outer corrections explicitly unless a future outer-state convergence metric is added.");
  params.setDocString("num_piso_iterations",
                      "The maximum number of additional inner pressure-correction-only PISO "
                      "stages performed without rebuilding the momentum matrix on each outer "
                      "correction. By default the reduced-pressure executioner performs this many "
                      "additional stages explicitly; early exit only occurs when explicit PISO "
                      "termination tolerances are provided.");
  params.addClassDescription(
      "PIMPLE solve object with explicit hooks for reduced-pressure sharp-interface face-flux "
      "predictors.");
  params.setDocString(
      "active_scalar_systems",
      "The solver system for each sharp-interface volume-fraction transport equation.");
  params.setDocString("should_solve_active_scalars",
                      "Whether to solve the volume-fraction transport equation(s).");
  params.setDocString("active_scalar_equation_relaxation",
                      "The relaxation used for the volume-fraction transport equation(s).");
  params.setDocString(
      "active_scalar_absolute_tolerance",
      "The absolute tolerance(s) on the normalized residual(s) of the volume-fraction "
      "equation(s).");
  params.addRangeCheckedParam<unsigned int>(
      "volume_fraction_subcycles",
      1,
      "volume_fraction_subcycles>0",
      "Number of full alpha transport subcycles used by the bounded volume-fraction update.");
  params.addRangeCheckedParam<Real>(
      "volume_fraction_max_courant",
      1.0,
      "volume_fraction_max_courant>0",
      "Maximum allowed alpha Courant number during subcycling. The executioner increases the "
      "alpha subcycle count as needed so the current transported volumetric flux satisfies this "
      "limit.");
  MooseEnum startup_pressure_initialization("none projection-only", "projection-only");
  params.addParam<MooseEnum>(
      "startup_pressure_initialization",
      startup_pressure_initialization,
      "Startup reduced-pressure initialization policy on the first time step. Use "
      "'projection-only' to apply startup continuity projection without overwriting the "
      "user-supplied reduced-pressure field, or 'none' to skip startup pressure "
      "cleanup entirely.");
  params.addRangeCheckedParam<unsigned int>(
      "startup_flux_corrections",
      1,
      "startup_flux_corrections>0",
      "Number of pressure-only startup cleanup / projection corrections applied when "
      "startup_pressure_initialization is not 'none'.");
  params.addParam<SolverSystemName>(
      "low_mach_mass_system",
      "",
      "Optional temporary-density system for the implicit conservative low-Mach mass update.");
  params.addParam<UserObjectName>(
      "low_mach_mass_flux",
      "",
      "LowMachImplicitMassFlux object that publishes the exact converged mass-equation face flux.");
  params.addParam<MultiMooseEnum>("low_mach_mass_petsc_options",
                                  Moose::PetscSupport::getCommonPetscFlags(),
                                  "Singleton PETSc options for the temporary-density equation.");
  params.addParam<MultiMooseEnum>(
      "low_mach_mass_petsc_options_iname",
      Moose::PetscSupport::getCommonPetscKeys(),
      "Names of PETSc name/value pairs for the temporary-density equation.");
  params.addParam<std::vector<std::string>>(
      "low_mach_mass_petsc_options_value",
      "Values of PETSc name/value pairs for the temporary-density equation.");
  params.addRangeCheckedParam<Real>(
      "low_mach_mass_l_tol",
      1e-10,
      "0.0<=low_mach_mass_l_tol & low_mach_mass_l_tol<1.0",
      "Relative linear-solver tolerance for each implicit mass correction.");
  params.addRangeCheckedParam<Real>(
      "low_mach_mass_l_abs_tol",
      1e-50,
      "0.0<low_mach_mass_l_abs_tol",
      "Absolute linear-solver tolerance, scaled by the mass-equation normalization factor.");
  params.addRangeCheckedParam<unsigned int>(
      "low_mach_mass_l_max_its",
      10000,
      "low_mach_mass_l_max_its>0",
      "Maximum linear iterations for each implicit mass correction.");
  params.addRangeCheckedParam<Real>(
      "low_mach_mass_absolute_tolerance",
      1e-8,
      "low_mach_mass_absolute_tolerance>0",
      "Outer-iteration convergence tolerance for the implicit mass update.");
  params.addRangeCheckedParam<unsigned int>(
      "low_mach_minimum_outer_iterations",
      2,
      "low_mach_minimum_outer_iterations>0",
      "Minimum number of complete low-Mach fixed-point sweeps. The paper uses two sweeps so the "
      "second material-indicator update consumes the velocity divergence produced by the first.");
  params.addRangeCheckedParam<Real>(
      "low_mach_fixed_point_tolerance",
      1e-6,
      "low_mach_fixed_point_tolerance>0",
      "Relative outer-state update tolerance for the coupled material-indicator, temporary-"
      "density, enthalpy, and pressure-corrected volumetric-flux fixed point.");
  params.addRangeCheckedParam<Real>(
      "low_mach_face_flux_absolute_tolerance",
      1e-12,
      "low_mach_face_flux_absolute_tolerance>0",
      "Absolute RMS tolerance for the pressure-corrected volumetric face-flux outer update. It is "
      "combined with low_mach_fixed_point_tolerance so nearly quiescent flow is not judged using "
      "a relative-to-zero flux norm.");
  params.addRangeCheckedParam<Real>(
      "low_mach_level_set_aitken_min_relaxation",
      0.05,
      "0<low_mach_level_set_aitken_min_relaxation & "
      "low_mach_level_set_aitken_min_relaxation<=1",
      "Lower bound for adaptive Aitken relaxation of the projected level-set outer update.");
  params.addRangeCheckedParam<Real>(
      "low_mach_level_set_aitken_max_relaxation",
      1.0,
      "0<low_mach_level_set_aitken_max_relaxation & "
      "low_mach_level_set_aitken_max_relaxation<=1",
      "Upper bound for adaptive Aitken relaxation of the projected level-set outer update.");
  params.addRangeCheckedParam<Real>(
      "low_mach_level_set_aitken_max_growth",
      1.5,
      "low_mach_level_set_aitken_max_growth>=1",
      "Maximum factor by which the projected level-set Aitken relaxation may increase between "
      "successive outer iterations.");
  params.addRangeCheckedParam<Real>(
      "low_mach_level_set_aitken_projection_amplification",
      2.0,
      "low_mach_level_set_aitken_projection_amplification>=1",
      "Maximum accepted ratio between the projected level-set update and the proposed relaxed "
      "update before the Aitken factor is reduced.");
  params.addParam<UserObjectName>(
      "low_mach_enthalpy_state",
      "",
      "LowMachEnthalpyNewtonState object that restores exact thermodynamic consistency after "
      "each linearized temperature solve.");
  params.addParam<UserObjectName>(
      "low_mach_divergence_source",
      "",
      "LowMachDivergenceSource object that constructs the pressure-constraint source from the "
      "converged enthalpy balance.");
  params.addParam<MooseFunctorName>(
      "low_mach_solid_drag",
      "low_mach_solid_drag",
      "Carman-Kozeny drag coefficient used implicitly by every momentum component.");
  params.addRangeCheckedParam<unsigned int>(
      "low_mach_enthalpy_max_iterations",
      5,
      "low_mach_enthalpy_max_iterations>0",
      "Maximum Newton/Picard corrections for the low-Mach enthalpy equation.");
  params.addRangeCheckedParam<Real>(
      "low_mach_enthalpy_tolerance",
      1e-8,
      "low_mach_enthalpy_tolerance>0",
      "Convergence tolerance for the normalized liquid-fraction L2 update.");
  params.addRangeCheckedParam<Real>(
      "low_mach_enthalpy_residual_tolerance",
      1e-8,
      "low_mach_enthalpy_residual_tolerance>0",
      "Convergence tolerance for the normalized nonlinear conservative-enthalpy residual.");
  params.addRangeCheckedParam<Real>(
      "low_mach_enthalpy_l_tol",
      1e-12,
      "0.0<low_mach_enthalpy_l_tol & low_mach_enthalpy_l_tol<1.0",
      "Relative linear-solver tolerance for each Newton temperature equation. This must be "
      "tighter than the nonlinear conservative-residual tolerance.");
  params.addParamNamesToGroup(
      "volume_fraction_subcycles volume_fraction_max_courant startup_pressure_initialization "
      "startup_flux_corrections num_pressure_nonorthogonal_correctors",
      "Volume Fraction Equations");
  params.addParamNamesToGroup("low_mach_mass_system low_mach_mass_flux low_mach_mass_petsc_options "
                              "low_mach_mass_petsc_options_iname low_mach_mass_petsc_options_value "
                              "low_mach_mass_l_tol low_mach_mass_l_abs_tol low_mach_mass_l_max_its "
                              "low_mach_mass_absolute_tolerance low_mach_minimum_outer_iterations "
                              "low_mach_fixed_point_tolerance "
                              "low_mach_face_flux_absolute_tolerance "
                              "low_mach_level_set_aitken_min_relaxation "
                              "low_mach_level_set_aitken_max_relaxation "
                              "low_mach_level_set_aitken_max_growth "
                              "low_mach_level_set_aitken_projection_amplification "
                              "low_mach_enthalpy_state low_mach_divergence_source "
                              "low_mach_solid_drag "
                              "low_mach_enthalpy_max_iterations low_mach_enthalpy_tolerance "
                              "low_mach_enthalpy_residual_tolerance low_mach_enthalpy_l_tol",
                              "Low-Mach Mass Equation");
  return params;
}

ReducedPressurePIMPLESolve::ReducedPressurePIMPLESolve(Executioner & ex)
  : PIMPLESolve(ex),
    _volume_fraction_subcycles(getParam<unsigned int>("volume_fraction_subcycles")),
    _volume_fraction_max_courant(getParam<Real>("volume_fraction_max_courant")),
    _startup_flux_corrections(getParam<unsigned int>("startup_flux_corrections")),
    _low_mach_mass_system_name(getParam<SolverSystemName>("low_mach_mass_system")),
    _has_low_mach_mass_system(!_low_mach_mass_system_name.empty()),
    _low_mach_mass_system_number(_has_low_mach_mass_system
                                     ? _problem.linearSysNum(_low_mach_mass_system_name)
                                     : libMesh::invalid_uint),
    _low_mach_mass_system(_has_low_mach_mass_system
                              ? &_problem.getLinearSystem(_low_mach_mass_system_number)
                              : nullptr),
    _low_mach_mass_flux(nullptr),
    _low_mach_mass_l_abs_tol(getParam<Real>("low_mach_mass_l_abs_tol")),
    _low_mach_mass_absolute_tolerance(getParam<Real>("low_mach_mass_absolute_tolerance")),
    _low_mach_minimum_outer_iterations(getParam<unsigned int>("low_mach_minimum_outer_iterations")),
    _low_mach_fixed_point_tolerance(getParam<Real>("low_mach_fixed_point_tolerance")),
    _low_mach_face_flux_absolute_tolerance(getParam<Real>("low_mach_face_flux_absolute_tolerance")),
    _low_mach_level_set_aitken_min_relaxation(
        getParam<Real>("low_mach_level_set_aitken_min_relaxation")),
    _low_mach_level_set_aitken_max_relaxation(
        getParam<Real>("low_mach_level_set_aitken_max_relaxation")),
    _low_mach_level_set_aitken_max_growth(getParam<Real>("low_mach_level_set_aitken_max_growth")),
    _low_mach_level_set_aitken_projection_amplification(
        getParam<Real>("low_mach_level_set_aitken_projection_amplification")),
    _low_mach_enthalpy_state(nullptr),
    _low_mach_divergence_source(nullptr),
    _low_mach_enthalpy_max_iterations(getParam<unsigned int>("low_mach_enthalpy_max_iterations")),
    _low_mach_enthalpy_tolerance(getParam<Real>("low_mach_enthalpy_tolerance")),
    _low_mach_enthalpy_residual_tolerance(getParam<Real>("low_mach_enthalpy_residual_tolerance")),
    _low_mach_enthalpy_l_tol(getParam<Real>("low_mach_enthalpy_l_tol"))
{
  _startup_pressure_initialization =
      getParam<MooseEnum>("startup_pressure_initialization").operator std::string();

  if (_low_mach_level_set_aitken_min_relaxation > _low_mach_level_set_aitken_max_relaxation)
    paramError("low_mach_level_set_aitken_min_relaxation",
               "The minimum Aitken relaxation must not exceed the maximum relaxation.");

  if (_pin_pressure)
    paramError("pin_pressure",
               "ReducedPressurePIMPLE supports pressure boundary-condition constraints only; "
               "pressure pinning is not supported.");

  const bool has_mass_flux = !getParam<UserObjectName>("low_mach_mass_flux").empty();
  if (_has_low_mach_mass_system != has_mass_flux)
    paramError(_has_low_mach_mass_system ? "low_mach_mass_flux" : "low_mach_mass_system",
               "'low_mach_mass_system' and 'low_mach_mass_flux' must be supplied together.");

  const bool has_enthalpy_state = !getParam<UserObjectName>("low_mach_enthalpy_state").empty();
  if (_has_low_mach_mass_system != has_enthalpy_state)
    paramError(_has_low_mach_mass_system ? "low_mach_enthalpy_state" : "low_mach_mass_system",
               "'low_mach_enthalpy_state' is required exactly when the low-Mach mass system is "
               "enabled.");

  const bool has_divergence_source =
      !getParam<UserObjectName>("low_mach_divergence_source").empty();
  if (_has_low_mach_mass_system != has_divergence_source)
    paramError(_has_low_mach_mass_system ? "low_mach_divergence_source" : "low_mach_mass_system",
               "'low_mach_divergence_source' is required exactly when the low-Mach mass system "
               "is enabled.");

  if (_has_low_mach_mass_system)
  {
    if (_num_iterations < _low_mach_minimum_outer_iterations)
      paramError("num_iterations",
                 "The low-Mach path requires at least ",
                 _low_mach_minimum_outer_iterations,
                 " outer fixed-point sweeps, but num_iterations is ",
                 _num_iterations,
                 ".");
    if (!_has_energy_system || !_should_solve_energy)
      paramError("low_mach_enthalpy_state",
                 "The low-Mach mass/enthalpy path requires an enabled fluid energy system.");
    if (!getParam<UserObjectName>("liquid_fraction_corrector").empty())
      paramError("liquid_fraction_corrector",
                 "Low-Mach enthalpy uses LowMachEnthalpyNewtonState; the laser-foam liquid "
                 "fraction corrector cannot be enabled simultaneously.");
    if (_has_solid_energy_system && _should_solve_solid_energy)
      paramError("should_solve_solid_energy",
                 "Conjugate solid-energy coupling has not been incorporated into the pre-flow "
                 "low-Mach enthalpy Newton loop.");
    if (_has_pm_radiation_systems && _should_solve_pm_radiation)
      paramError("should_solve_pm_radiation",
                 "Participating-media radiation has not been incorporated into the pre-flow "
                 "low-Mach enthalpy Newton loop.");
    if (_cht.enabled())
      paramError("low_mach_enthalpy_state",
                 "Conjugate heat transfer has not been incorporated into the pre-flow low-Mach "
                 "enthalpy Newton loop.");

    _systems_to_solve.push_back(_low_mach_mass_system);
    _low_mach_mass_system->system().prefix_with_name(false);

    const auto & petsc_flags = getParam<MultiMooseEnum>("low_mach_mass_petsc_options");
    const auto & petsc_pairs = getParam<MooseEnumItem, std::string>(
        "low_mach_mass_petsc_options_iname", "low_mach_mass_petsc_options_value");
    Moose::PetscSupport::addPetscFlagsToPetscOptions(
        petsc_flags, "", *this, _low_mach_mass_petsc_options);
    Moose::PetscSupport::addPetscPairsToPetscOptions(
        petsc_pairs, _problem.mesh().dimension(), "", *this, _low_mach_mass_petsc_options);

    _low_mach_mass_linear_control.real_valued_data["rel_tol"] =
        getParam<Real>("low_mach_mass_l_tol");
    _low_mach_mass_linear_control.real_valued_data["abs_tol"] = _low_mach_mass_l_abs_tol;
    _low_mach_mass_linear_control.int_valued_data["max_its"] =
        getParam<unsigned int>("low_mach_mass_l_max_its");
    _energy_linear_control.real_valued_data["rel_tol"] = _low_mach_enthalpy_l_tol;
  }
}

void
ReducedPressurePIMPLESolve::initialSetup()
{
  PIMPLESolve::initialSetup();

  if (_has_low_mach_mass_system)
  {
    _low_mach_mass_flux = const_cast<LowMachImplicitMassFlux *>(
        &getUserObject<LowMachImplicitMassFlux>("low_mach_mass_flux"));
    if (_low_mach_mass_flux->systemName() != _low_mach_mass_system_name)
      paramError("low_mach_mass_flux",
                 "The LowMachImplicitMassFlux object is linked to system '",
                 _low_mach_mass_flux->systemName(),
                 "', but the executioner is configured to solve '",
                 _low_mach_mass_system_name,
                 "'.");

    _low_mach_enthalpy_state = const_cast<LowMachEnthalpyNewtonState *>(
        &getUserObject<LowMachEnthalpyNewtonState>("low_mach_enthalpy_state"));
    if (_low_mach_enthalpy_state->systemName() != _energy_system->name())
      paramError("low_mach_enthalpy_state",
                 "The LowMachEnthalpyNewtonState object is linked to system '",
                 _low_mach_enthalpy_state->systemName(),
                 "', but the configured energy system is '",
                 _energy_system->name(),
                 "'.");

    _low_mach_divergence_source = const_cast<LowMachDivergenceSource *>(
        &getUserObject<LowMachDivergenceSource>("low_mach_divergence_source"));
    if (_low_mach_divergence_source->systemName() != _energy_system->name())
      paramError("low_mach_divergence_source",
                 "The LowMachDivergenceSource object is linked to system '",
                 _low_mach_divergence_source->systemName(),
                 "', but the configured energy system is '",
                 _energy_system->name(),
                 "'.");

    validateLowMachEnthalpyOperators();
  }
}

bool
ReducedPressurePIMPLESolve::startupPressureInitializationEnabled() const
{
  return _startup_pressure_initialization != "none";
}

bool
ReducedPressurePIMPLESolve::shouldRunStartupInitialization() const
{
  return startupPressureInitializationEnabled() && _problem.timeStep() == 1;
}

bool
ReducedPressurePIMPLESolve::solvesVolumeFraction() const
{
  return _has_active_scalar_systems && _should_solve_active_scalars;
}

ConservativeSharpInterfaceVOFMULESCorrector *
ReducedPressurePIMPLESolve::sharpInterfaceVOFCorrector(const SolverSystemName & system_name) const
{
  std::vector<UserObject *> objs;
  _problem.theWarehouse()
      .query()
      .condition<AttribSystem>("UserObject")
      .condition<AttribThread>(0)
      .queryInto(objs);

  ConservativeSharpInterfaceVOFMULESCorrector * corrector_match = nullptr;
  for (const auto & obj : objs)
    if (auto * corrector = dynamic_cast<ConservativeSharpInterfaceVOFMULESCorrector *>(obj);
        corrector && corrector->systemName() == system_name)
    {
      if (corrector_match)
        mooseError(
            "ReducedPressurePIMPLESolve found multiple ConservativeSharpInterfaceVOFMULESCorrector "
            "objects for system '",
            system_name,
            "'.");
      corrector_match = corrector;
    }

  return corrector_match;
}

MassConservedLevelSetCorrector *
ReducedPressurePIMPLESolve::massConservedLevelSetCorrector(
    const SolverSystemName & system_name) const
{
  std::vector<UserObject *> objs;
  _problem.theWarehouse()
      .query()
      .condition<AttribSystem>("UserObject")
      .condition<AttribThread>(0)
      .queryInto(objs);

  MassConservedLevelSetCorrector * corrector_match = nullptr;
  for (const auto & obj : objs)
    if (auto * corrector = dynamic_cast<MassConservedLevelSetCorrector *>(obj);
        corrector && corrector->systemName() == system_name)
    {
      if (corrector_match)
        mooseError("ReducedPressurePIMPLESolve found multiple MassConservedLevelSetCorrector "
                   "objects for system '",
                   system_name,
                   "'.");
      corrector_match = corrector;
    }

  return corrector_match;
}

void
ReducedPressurePIMPLESolve::preSolveSetup(const SolverParams & /* solver_params */)
{
  if (auto * sharp_rc = sharpInterfaceRC())
    sharp_rc->setSuppressExplicitHydrostaticPressureFlux(false);
}

void
ReducedPressurePIMPLESolve::addIterationResiduals(ResidualStorage & residual_storage)
{
  PIMPLESolve::addIterationResiduals(residual_storage);

  if (_has_low_mach_mass_system)
  {
    _low_mach_mass_residual_index = residual_storage.ns_residuals.size();
    residual_storage.ns_residuals.emplace_back(0, std::numeric_limits<Real>::max());
    residual_storage.ns_abs_tols.push_back(_low_mach_mass_absolute_tolerance);
    residual_storage.converged = false;

    _low_mach_liquid_fraction_residual_index = residual_storage.ns_residuals.size();
    residual_storage.ns_residuals.emplace_back(0, std::numeric_limits<Real>::max());
    residual_storage.ns_abs_tols.push_back(_low_mach_enthalpy_tolerance);

    _low_mach_fixed_point_residual_index = residual_storage.ns_residuals.size();
    residual_storage.ns_residuals.emplace_back(0, 1.0);
    residual_storage.ns_abs_tols.push_back(_low_mach_fixed_point_tolerance);
  }
}

void
ReducedPressurePIMPLESolve::initializeSolveLoop(const SolverParams & solver_params)
{
  resetVOFTransportStateForNewSolve();

  if (_has_low_mach_mass_system && _problem.timeStep() == 1)
  {
    _problem.execute(EXEC_NONLINEAR);
    _low_mach_mass_flux->synchronizeDensityToEOS();
    synchronizeSystemState(*_low_mach_mass_system);
  }

  if (shouldRunStartupInitialization())
    initializeConsistentStartupState();

  if (_should_solve_pressure)
  {
    initializeStartupPressureField(solver_params);
    commitAcceptedVOFTransportHistoryIfNeeded();
  }

  if (_has_low_mach_mass_system)
    captureLowMachOuterState();
}

void
ReducedPressurePIMPLESolve::preMomentumPressureIteration(ResidualStorage & residual_storage,
                                                         const SolverParams & solver_params)
{
  if (_has_low_mach_mass_system)
  {
    const Real outer_update = _current_outer_iteration == 1 ? 1.0 : updateLowMachOuterState();
    residual_storage.ns_residuals[_low_mach_fixed_point_residual_index] = std::make_pair(
        0,
        _current_outer_iteration < _low_mach_minimum_outer_iterations ? std::max(1.0, outer_update)
                                                                      : outer_update);
  }

  advanceOuterIterationHistories();

  if (solvesVolumeFraction())
    solveVolumeFractionBeforeFlowCorrection(residual_storage, solver_params);

  if (_has_low_mach_mass_system)
  {
    solveImplicitLowMachMassEquation(residual_storage, solver_params);
    solveImplicitLowMachEnthalpyEquation(residual_storage, solver_params);
    _low_mach_divergence_source->publish();
  }
}

void
ReducedPressurePIMPLESolve::solveImplicitLowMachMassEquation(ResidualStorage & residual_storage,
                                                             const SolverParams & solver_params)
{
  mooseAssert(_low_mach_mass_system, "The low-Mach mass system must be available.");
  mooseAssert(_low_mach_mass_flux, "The low-Mach mass-flux publisher must be linked.");

  _problem.execute(EXEC_NONLINEAR);
  Moose::PetscSupport::petscSetOptions(_low_mach_mass_petsc_options, solver_params);

  auto & linear_system = libMesh::cast_ref<LinearImplicitSystem &>(_low_mach_mass_system->system());
  auto & solution = *linear_system.solution;
  auto previous_solution = solution.clone();
  const auto mass_residual = solveAdvectedSystem(_low_mach_mass_system_number,
                                                 *_low_mach_mass_system,
                                                 1.0,
                                                 _low_mach_mass_linear_control,
                                                 _low_mach_mass_l_abs_tol,
                                                 1.0,
                                                 std::numeric_limits<Real>::min());
  const Real solution_norm = solution.l2_norm();
  const Real relative_update =
      solution.l2_norm_diff(*previous_solution) /
      std::max(solution_norm, std::sqrt(std::numeric_limits<Real>::epsilon()));

  _console << " Low-Mach implicit upwind update: " << COLOR_GREEN << relative_update
           << COLOR_DEFAULT << std::endl;

  _low_mach_mass_flux->publishConvergedMassFlux();
  _problem.execute(EXEC_NONLINEAR);

  residual_storage.ns_residuals[_low_mach_mass_residual_index] =
      std::make_pair(mass_residual.first, relative_update);
}

void
ReducedPressurePIMPLESolve::solveImplicitLowMachEnthalpyEquation(ResidualStorage & residual_storage,
                                                                 const SolverParams & solver_params)
{
  mooseAssert(_energy_system, "The low-Mach temperature system must be available.");
  mooseAssert(_low_mach_enthalpy_state, "The low-Mach enthalpy state corrector must be linked.");

  Moose::PetscSupport::petscSetOptions(_energy_petsc_options, solver_params);
  unsigned int total_linear_iterations = 0;
  unsigned int nonlinear_iterations = 0;
  Real liquid_fraction_update = std::numeric_limits<Real>::max();
  _low_mach_divergence_source->beginFixedPointEnthalpySolve();
  auto nonlinear_residual = assembleLowMachEnthalpyResidual();
  _low_mach_divergence_source->finishFixedPointEnthalpyInitialization();
  std::pair<unsigned int, Real> energy_residual{0, std::numeric_limits<Real>::max()};

  for (const auto nonlinear_iteration : make_range(_low_mach_enthalpy_max_iterations))
  {
    if (nonlinear_residual.converged(_low_mach_enthalpy_residual_tolerance))
    {
      liquid_fraction_update = 0.0;
      break;
    }

    nonlinear_iterations = nonlinear_iteration + 1;
    _problem.execute(EXEC_NONLINEAR);
    _low_mach_enthalpy_state->beginNewtonIteration();
    _low_mach_divergence_source->beginEnthalpyAssembly();

    energy_residual = solveAdvectedSystem(_energy_sys_number,
                                          *_energy_system,
                                          1.0,
                                          _energy_linear_control,
                                          _energy_l_abs_tol,
                                          1.0,
                                          std::numeric_limits<Real>::lowest());
    total_linear_iterations += energy_residual.first;

    _low_mach_enthalpy_state->captureNewtonUpdate();

    constexpr Real reduction = 0.5;
    constexpr Real sufficient_decrease = 1e-4;
    constexpr unsigned int maximum_backtracks = 12;
    Real step_length = 1.0;
    LowMachEnthalpyResidual trial_residual{std::numeric_limits<Real>::max(), 1.0, 1.0};
    unsigned int backtracks = 0;
    while (true)
    {
      liquid_fraction_update = _low_mach_enthalpy_state->restoreExactState(step_length);
      _problem.execute(EXEC_NONLINEAR);
      trial_residual = assembleLowMachEnthalpyResidual();

      if (trial_residual.converged(_low_mach_enthalpy_residual_tolerance) ||
          trial_residual.absolute <=
              (1.0 - sufficient_decrease * step_length) * nonlinear_residual.absolute)
        break;
      if (backtracks == maximum_backtracks)
        mooseError("The low-Mach enthalpy line search failed to reduce the absolute nonlinear "
                   "residual from ",
                   nonlinear_residual.absolute,
                   " (normalized: ",
                   nonlinear_residual.normalized(),
                   "). Final trial residual: ",
                   trial_residual.absolute,
                   " (normalized: ",
                   trial_residual.normalized(),
                   ")",
                   " at step length ",
                   step_length,
                   ".");

      step_length *= reduction;
      ++backtracks;
    }
    nonlinear_residual = trial_residual;

    _console << " Low-Mach enthalpy Newton correction " << nonlinear_iteration + 1 << ": "
             << COLOR_GREEN << liquid_fraction_update << COLOR_DEFAULT
             << ", nonlinear residual: " << nonlinear_residual.normalized()
             << ", step length: " << step_length << std::endl;

    if (liquid_fraction_update <= _low_mach_enthalpy_tolerance &&
        nonlinear_residual.converged(_low_mach_enthalpy_residual_tolerance))
      break;
  }

  if (liquid_fraction_update > _low_mach_enthalpy_tolerance ||
      !nonlinear_residual.converged(_low_mach_enthalpy_residual_tolerance))
    mooseError("The low-Mach enthalpy solve did not converge in ",
               _low_mach_enthalpy_max_iterations,
               " globalized Newton corrections. Final normalized liquid-fraction update: ",
               liquid_fraction_update,
               "; final normalized nonlinear residual: ",
               nonlinear_residual.normalized(),
               ".");

  residual_storage.ns_residuals[residual_storage.energy_index] =
      std::make_pair(total_linear_iterations,
                     nonlinear_residual.converged(_low_mach_enthalpy_residual_tolerance)
                         ? 0.0
                         : nonlinear_residual.normalized());
  residual_storage.ns_residuals[_low_mach_liquid_fraction_residual_index] =
      std::make_pair(nonlinear_iterations, liquid_fraction_update);
}

ReducedPressurePIMPLESolve::LowMachEnthalpyResidual
ReducedPressurePIMPLESolve::assembleLowMachEnthalpyResidual()
{
  mooseAssert(_energy_system, "The low-Mach temperature system must be available.");
  mooseAssert(_low_mach_divergence_source, "The low-Mach divergence source must be available.");

  _problem.setCurrentLinearSystem(_energy_sys_number);
  auto & system = libMesh::cast_ref<LinearImplicitSystem &>(_energy_system->system());
  auto & solution = *system.solution;
  auto & matrix = *system.matrix;
  auto & rhs = *system.rhs;

  _low_mach_divergence_source->beginEnthalpyAssembly();
  _problem.computeLinearSystemSys(system, matrix, rhs, true);

  auto residual = solution.zero_clone();
  matrix.vector_mult(*residual, solution);
  const Real algebraic_scale = residual->l2_norm() + rhs.l2_norm();
  residual->add(-1.0, rhs);

  const Real normalization = NS::FV::computeNormalizationFactor(solution, matrix, rhs);
  mooseAssert(normalization > 0.0, "The enthalpy residual normalization must be positive.");
  return {residual->l2_norm(), normalization, algebraic_scale};
}

void
ReducedPressurePIMPLESolve::validateLowMachEnthalpyOperators() const
{
  if (_low_mach_enthalpy_state->materialFractionName() !=
      _low_mach_divergence_source->materialFractionName())
    mooseError("Low-Mach enthalpy state '",
               _low_mach_enthalpy_state->name(),
               "' and divergence source '",
               _low_mach_divergence_source->name(),
               "' must consume the same gas/PCM material-fraction functor.");

  if (solvesVolumeFraction())
  {
    ConservativeSharpInterfaceVOFMULESCorrector * material_indicator_corrector = nullptr;
    MassConservedLevelSetCorrector * level_set_corrector = nullptr;
    for (const auto * active_scalar_system : _active_scalar_systems)
    {
      if (auto * corrector = sharpInterfaceVOFCorrector(active_scalar_system->name()))
      {
        if (material_indicator_corrector || level_set_corrector)
          mooseError("The low-Mach enthalpy path requires exactly one transported gas/PCM "
                     "material indicator, but multiple interface correctors are active.");
        material_indicator_corrector = corrector;
      }
      if (auto * corrector = massConservedLevelSetCorrector(active_scalar_system->name()))
      {
        if (material_indicator_corrector || level_set_corrector)
          mooseError("The low-Mach enthalpy path requires exactly one transported gas/PCM "
                     "material indicator, but multiple interface correctors are active.");
        level_set_corrector = corrector;
      }
    }

    if (!material_indicator_corrector && !level_set_corrector)
      mooseError("The low-Mach enthalpy path has active scalar transport enabled, but no "
                 "VOF or mass-conserved level-set corrector transports its gas/PCM material "
                 "indicator.");

    if (material_indicator_corrector)
    {
      if (material_indicator_corrector->variableName() !=
          _low_mach_enthalpy_state->materialFractionName())
        mooseError("Transported material indicator '",
                   material_indicator_corrector->variableName(),
                   "' must be the material-fraction functor used by the low-Mach thermodynamic "
                   "state, which is currently '",
                   _low_mach_enthalpy_state->materialFractionName(),
                   "'.");

      if (material_indicator_corrector->sourceSpName() != _low_mach_divergence_source->sourceName())
        mooseError("The material-indicator equation must use the low-Mach divergence source '",
                   _low_mach_divergence_source->sourceName(),
                   "' as its implicit source coefficient so it discretizes "
                   "dH/dt + div(H u) - H div(u) = 0.");
    }
    else if (level_set_corrector->materialFractionName() !=
             _low_mach_enthalpy_state->materialFractionName())
      mooseError("Mass-conserved level set '",
                 level_set_corrector->variableName(),
                 "' reconstructs material-fraction functor '",
                 level_set_corrector->materialFractionName(),
                 "', but the low-Mach thermodynamic state uses '",
                 _low_mach_enthalpy_state->materialFractionName(),
                 "'.");

    std::vector<LinearFVFluxKernel *> mass_flux_kernels;
    _problem.theWarehouse()
        .query()
        .condition<AttribThread>(0)
        .condition<AttribSysNum>(_low_mach_mass_system_number)
        .condition<AttribSystem>("LinearFVFluxKernel")
        .queryInto(mass_flux_kernels);
    const LinearFVLowMachMassAdvection * mass_advection = nullptr;
    for (const auto * kernel : mass_flux_kernels)
      if (const auto * candidate = dynamic_cast<const LinearFVLowMachMassAdvection *>(kernel))
      {
        if (mass_advection)
          mooseError("The low-Mach temporary-density system must contain exactly one "
                     "LinearFVLowMachMassAdvection kernel.");
        mass_advection = candidate;
      }
    if (!mass_advection)
      mooseError("The low-Mach temporary-density system contains no "
                 "LinearFVLowMachMassAdvection kernel.");

    const MooseFunctorName & indicator_flux = material_indicator_corrector
                                                  ? material_indicator_corrector->faceFluxName()
                                                  : level_set_corrector->volumetricFaceFluxName();
    if (indicator_flux != mass_advection->volumetricFaceFluxName())
      mooseError("The material-indicator and temporary-density equations must use the same "
                 "volumetric face flux. They currently use '",
                 indicator_flux,
                 "' and '",
                 mass_advection->volumetricFaceFluxName(),
                 "', respectively.");
  }

  std::vector<LinearFVFluxKernel *> flux_kernels;
  _problem.theWarehouse()
      .query()
      .condition<AttribThread>(0)
      .condition<AttribSysNum>(_energy_sys_number)
      .condition<AttribSystem>("LinearFVFluxKernel")
      .queryInto(flux_kernels);

  unsigned int enthalpy_advection_kernels = 0;
  for (const auto * kernel : flux_kernels)
    if (const auto * enthalpy_advection =
            dynamic_cast<const LinearFVLowMachEnthalpyAdvection *>(kernel))
    {
      ++enthalpy_advection_kernels;
      if (!enthalpy_advection->usesMassFlux(_low_mach_mass_flux->massFluxName()))
        mooseError("Low-Mach enthalpy advection kernel '",
                   enthalpy_advection->name(),
                   "' must consume the non-integrated mass-flux functor published by '",
                   _low_mach_mass_flux->name(),
                   "'.");
    }

  if (!enthalpy_advection_kernels)
    mooseError("The configured low-Mach energy system contains no "
               "LinearFVLowMachEnthalpyAdvection kernel.");

  for (const auto * momentum_system : _momentum_systems)
  {
    std::vector<LinearFVFluxKernel *> momentum_kernels;
    _problem.theWarehouse()
        .query()
        .condition<AttribThread>(0)
        .condition<AttribSysNum>(momentum_system->number())
        .condition<AttribSystem>("LinearFVFluxKernel")
        .queryInto(momentum_kernels);

    unsigned int momentum_flux_kernels = 0;
    for (const auto * kernel : momentum_kernels)
      if (const auto * momentum_flux = dynamic_cast<const LinearWCNSFVMomentumFlux *>(kernel))
      {
        ++momentum_flux_kernels;
        if (!momentum_flux->usesMassFlux(_low_mach_mass_flux->massFluxName()))
          mooseError("Low-Mach momentum flux kernel '",
                     momentum_flux->name(),
                     "' must consume the non-integrated mass-flux functor published by '",
                     _low_mach_mass_flux->name(),
                     "'.");
      }

    if (!momentum_flux_kernels)
      mooseError("Low-Mach momentum system '",
                 momentum_system->name(),
                 "' contains no LinearWCNSFVMomentumFlux kernel.");

    std::vector<LinearFVElementalKernel *> momentum_elemental_kernels;
    _problem.theWarehouse()
        .query()
        .condition<AttribThread>(0)
        .condition<AttribSysNum>(momentum_system->number())
        .condition<AttribSystem>("LinearFVElementalKernel")
        .queryInto(momentum_elemental_kernels);

    unsigned int momentum_time_kernels = 0;
    unsigned int solid_drag_kernels = 0;
    for (const auto * kernel : momentum_elemental_kernels)
    {
      if (const auto * time_derivative = dynamic_cast<const LinearFVTimeDerivative *>(kernel))
      {
        ++momentum_time_kernels;
        if (!time_derivative->usesConservativeFactor(_low_mach_mass_flux->densityVariableName()))
          mooseError("Low-Mach momentum time kernel '",
                     time_derivative->name(),
                     "' must use the temporary-density functor with its old time state on the "
                     "right-hand side.");
      }
      if (const auto * solid_drag = dynamic_cast<const LinearFVLowMachSolidDrag *>(kernel))
      {
        ++solid_drag_kernels;
        if (!solid_drag->usesCoefficient(getParam<MooseFunctorName>("low_mach_solid_drag")))
          mooseError("Low-Mach solid-drag kernel '",
                     solid_drag->name(),
                     "' must use the configured Carman-Kozeny drag coefficient.");
      }
    }

    if (!momentum_time_kernels)
      mooseError("Low-Mach momentum system '",
                 momentum_system->name(),
                 "' contains no LinearFVTimeDerivative kernel.");
    if (solid_drag_kernels != 1)
      mooseError("Low-Mach momentum system '",
                 momentum_system->name(),
                 "' must contain exactly one LinearFVLowMachSolidDrag kernel.");
  }

  std::vector<LinearFVElementalKernel *> pressure_kernels;
  _problem.theWarehouse()
      .query()
      .condition<AttribThread>(0)
      .condition<AttribSysNum>(_pressure_sys_number)
      .condition<AttribSystem>("LinearFVElementalKernel")
      .queryInto(pressure_kernels);

  unsigned int divergence_kernels = 0;
  for (const auto * kernel : pressure_kernels)
    if (const auto * divergence = dynamic_cast<const LinearFVLowMachDivergence *>(kernel))
    {
      ++divergence_kernels;
      if (!divergence->usesSource(_low_mach_divergence_source->sourceName()))
        mooseError("Low-Mach divergence kernel '",
                   divergence->name(),
                   "' must consume the source functor published by '",
                   _low_mach_divergence_source->name(),
                   "'.");
    }

  if (!divergence_kernels)
    mooseError("The configured pressure system contains no LinearFVLowMachDivergence kernel.");
}

void
ReducedPressurePIMPLESolve::resetVOFTransportStateForNewSolve() const
{
  if (solvesVolumeFraction())
    if (auto * sharp_rc = sharpInterfaceRC())
      sharp_rc->clearVOFTransportState();
}

void
ReducedPressurePIMPLESolve::initializeConsistentStartupState()
{
  for (auto * system : _momentum_systems)
    synchronizeSystemState(*system);
  synchronizeSystemState(_pressure_system);
  for (auto * system : _active_scalar_systems)
    synchronizeSystemState(*system);
  if (_low_mach_mass_system)
    synchronizeSystemState(*_low_mach_mass_system);

  _problem.execute(EXEC_NONLINEAR);
}

void
ReducedPressurePIMPLESolve::captureLowMachOuterState()
{
  _low_mach_previous_outer_solutions.clear();
  const auto capture_solution = [this](LinearSystem & system)
  {
    auto & solution = system.solution();
    solution.close();
    _low_mach_previous_outer_solutions.push_back(solution.clone());
  };

  for (auto * system : _active_scalar_systems)
    capture_solution(*system);
  capture_solution(*_low_mach_mass_system);
  capture_solution(*_energy_system);

  _low_mach_previous_outer_face_flux.clear();
  if (_rc_uo)
    for (auto face = _problem.mesh().ownedFaceInfoBegin();
         face != _problem.mesh().ownedFaceInfoEnd();
         ++face)
      _low_mach_previous_outer_face_flux.emplace((*face)->id(),
                                                 _rc_uo->getVolumetricFaceFlux(**face));

  _level_set_previous_aitken_residuals.clear();
  _level_set_previous_aitken_residuals.resize(_active_scalar_systems.size());
  _level_set_aitken_relaxation.assign(_active_scalar_systems.size(),
                                      _low_mach_level_set_aitken_max_relaxation);
  _level_set_raw_fixed_point_update.assign(_active_scalar_systems.size(),
                                           std::numeric_limits<Real>::quiet_NaN());
}

Real
ReducedPressurePIMPLESolve::updateLowMachOuterState()
{
  std::vector<LinearSystem *> systems;
  systems.reserve(_active_scalar_systems.size() + 2);
  systems.insert(systems.end(), _active_scalar_systems.begin(), _active_scalar_systems.end());
  systems.push_back(_low_mach_mass_system);
  systems.push_back(_energy_system);
  mooseAssert(systems.size() == _low_mach_previous_outer_solutions.size(),
              "The low-Mach outer-state snapshot has the wrong number of systems.");

  Real maximum_update = 0.0;
  for (const auto i : index_range(systems))
  {
    auto & solution = systems[i]->solution();
    solution.close();
    const Real previous_norm = _low_mach_previous_outer_solutions[i]->l2_norm();
    const Real current_norm = solution.l2_norm();
    const Real absolute_update = solution.l2_norm_diff(*_low_mach_previous_outer_solutions[i]);
    const Real relative_update =
        absolute_update /
        std::max({current_norm, previous_norm, std::sqrt(std::numeric_limits<Real>::epsilon())});
    maximum_update = std::max(maximum_update, relative_update);
    _console << " Low-Mach outer " << systems[i]->name() << " update: relative=" << COLOR_GREEN
             << relative_update << COLOR_DEFAULT << ", absolute L2=" << absolute_update
             << std::endl;
    if (i < _active_scalar_systems.size() && std::isfinite(_level_set_raw_fixed_point_update[i]))
    {
      maximum_update = std::max(maximum_update, _level_set_raw_fixed_point_update[i]);
      _console << " Low-Mach outer " << systems[i]->name()
               << " raw projected fixed-point residual: " << COLOR_GREEN
               << _level_set_raw_fixed_point_update[i] << COLOR_DEFAULT << std::endl;
    }
    *_low_mach_previous_outer_solutions[i] = solution;
    _low_mach_previous_outer_solutions[i]->close();
  }

  if (_rc_uo)
  {
    Real local_flux_update_squared = 0.0;
    Real local_flux_squared = 0.0;
    Real local_previous_flux_squared = 0.0;
    dof_id_type local_face_count = 0;
    for (auto face = _problem.mesh().ownedFaceInfoBegin();
         face != _problem.mesh().ownedFaceInfoEnd();
         ++face)
    {
      const auto face_id = (*face)->id();
      const Real flux = _rc_uo->getVolumetricFaceFlux(**face);
      const auto previous_flux = _low_mach_previous_outer_face_flux.find(face_id);
      mooseAssert(previous_flux != _low_mach_previous_outer_face_flux.end(),
                  "A pressure-corrected face flux is missing from the low-Mach snapshot.");
      local_flux_update_squared += Utility::pow<2>(flux - previous_flux->second);
      local_flux_squared += Utility::pow<2>(flux);
      local_previous_flux_squared += Utility::pow<2>(previous_flux->second);
      ++local_face_count;
      previous_flux->second = flux;
    }
    _communicator.sum(local_flux_update_squared);
    _communicator.sum(local_flux_squared);
    _communicator.sum(local_previous_flux_squared);
    _communicator.sum(local_face_count);
    mooseAssert(local_face_count, "The low-Mach mesh contains no owned faces.");

    const Real inverse_face_count = 1.0 / local_face_count;
    const Real flux_update_rms = std::sqrt(local_flux_update_squared * inverse_face_count);
    const Real flux_rms = std::sqrt(local_flux_squared * inverse_face_count);
    const Real previous_flux_rms = std::sqrt(local_previous_flux_squared * inverse_face_count);
    // update <= relative_tolerance * reference + absolute_tolerance is equivalent to comparing
    // this combined normalized update to the relative fixed-point tolerance.
    const Real normalized_flux_update = flux_update_rms / (std::max(flux_rms, previous_flux_rms) +
                                                           _low_mach_face_flux_absolute_tolerance /
                                                               _low_mach_fixed_point_tolerance);
    maximum_update = std::max(maximum_update, normalized_flux_update);
    _console << " Low-Mach outer pressure-corrected face-flux update: combined=" << COLOR_GREEN
             << normalized_flux_update << COLOR_DEFAULT << ", absolute RMS=" << flux_update_rms
             << ", state RMS=" << flux_rms << std::endl;
  }

  _console << " Low-Mach outer fixed-point update: " << COLOR_GREEN << maximum_update
           << COLOR_DEFAULT << std::endl;
  return maximum_update;
}

void
ReducedPressurePIMPLESolve::applyLevelSetAitkenRelaxation(
    const unsigned int system_index,
    LinearSystem & system,
    MassConservedLevelSetCorrector & corrector)
{
  mooseAssert(system_index < _active_scalar_systems.size(),
              "The level-set system index is out of range.");
  mooseAssert(system_index < _low_mach_previous_outer_solutions.size(),
              "The previous level-set outer state is unavailable.");
  mooseAssert(system_index < _level_set_previous_aitken_residuals.size(),
              "The previous level-set Aitken residual is unavailable.");
  mooseAssert(system_index < _level_set_raw_fixed_point_update.size(),
              "The raw level-set fixed-point residual storage is unavailable.");

  auto & solution = system.solution();
  solution.close();
  const auto & previous_outer_solution = *_low_mach_previous_outer_solutions[system_index];
  auto residual = solution.clone();
  residual->add(-1.0, previous_outer_solution);
  residual->close();

  const Real solution_scale =
      std::max(solution.l2_norm(), std::sqrt(std::numeric_limits<Real>::epsilon()));
  const Real raw_update_norm = residual->l2_norm();
  const Real raw_update = raw_update_norm / solution_scale;
  _level_set_raw_fixed_point_update[system_index] = raw_update;
  Real relaxation = _level_set_aitken_relaxation[system_index];

  if (const auto & previous_residual = _level_set_previous_aitken_residuals[system_index])
  {
    auto residual_change = residual->clone();
    residual_change->add(-1.0, *previous_residual);
    residual_change->close();

    const Real denominator = residual_change->dot(*residual_change);
    const Real residual_scale =
        residual->dot(*residual) + previous_residual->dot(*previous_residual);
    if (denominator > std::numeric_limits<Real>::epsilon() *
                          std::max(residual_scale, std::numeric_limits<Real>::min()))
    {
      const Real aitken_relaxation =
          -relaxation * previous_residual->dot(*residual_change) / denominator;
      if (std::isfinite(aitken_relaxation))
      {
        const Real maximum_growth_relaxation =
            std::min(_low_mach_level_set_aitken_max_relaxation,
                     relaxation * _low_mach_level_set_aitken_max_growth);
        relaxation = std::clamp(aitken_relaxation,
                                _low_mach_level_set_aitken_min_relaxation,
                                maximum_growth_relaxation);
      }
    }
  }

  _level_set_previous_aitken_residuals[system_index] = residual->clone();

  const auto project_relaxed_trial =
      [&system, &solution, &previous_outer_solution, &residual, &corrector](const Real omega)
  {
    solution = previous_outer_solution;
    solution.add(omega, *residual);
    solution.close();

    auto & current_local_solution = *system.system().current_local_solution;
    solution.localize(current_local_solution, system.dofMap().get_send_list());
    system.setSolution(current_local_solution);
    system.computeGradients();

    // The transported candidate was already fully redistanced. Repeating redistancing after
    // blending can move the zero contour and undo the relaxed interface displacement. Restore only
    // the exact nonlinear mass constraint here; the next transport candidate is fully redistanced.
    corrector.correctMassOnly();
  };

  Real accepted_update_norm = raw_update_norm;
  bool projection_limited = false;
  bool used_projected_candidate = false;
  if (relaxation < 1.0)
  {
    while (true)
    {
      project_relaxed_trial(relaxation);
      solution.close();
      accepted_update_norm = solution.l2_norm_diff(previous_outer_solution);

      const Real permitted_projected_update =
          _low_mach_level_set_aitken_projection_amplification * relaxation * raw_update_norm;
      if (accepted_update_norm <= permitted_projected_update)
        break;

      projection_limited = true;
      if (relaxation <= _low_mach_level_set_aitken_min_relaxation)
      {
        // Prefer the fully projected fixed-point candidate to a nominally relaxed state that the
        // nonlinear projection moves farther than the proposed relaxed update.
        solution = previous_outer_solution;
        solution.add(1.0, *residual);
        solution.close();
        auto & current_local_solution = *system.system().current_local_solution;
        solution.localize(current_local_solution, system.dofMap().get_send_list());
        system.setSolution(current_local_solution);
        system.computeGradients();
        relaxation = 1.0;
        accepted_update_norm = raw_update_norm;
        used_projected_candidate = true;
        break;
      }

      relaxation = std::max(_low_mach_level_set_aitken_min_relaxation, 0.5 * relaxation);
    }
  }

  _level_set_aitken_relaxation[system_index] = relaxation;
  solution.close();
  const Real accepted_update =
      accepted_update_norm /
      std::max(solution.l2_norm(), std::sqrt(std::numeric_limits<Real>::epsilon()));
  _console << " Low-Mach level-set Aitken relaxation: omega=" << COLOR_GREEN << relaxation
           << COLOR_DEFAULT << ", raw update=" << raw_update
           << ", accepted update=" << accepted_update;
  if (projection_limited)
    _console << " (projection safeguard active)";
  if (used_projected_candidate)
    _console << " (using fully projected candidate)";
  _console << std::endl;
}

void
ReducedPressurePIMPLESolve::commitAcceptedVOFTransportHistoryIfNeeded() const
{
  if (solvesVolumeFraction())
    if (auto * sharp_rc = sharpInterfaceRC())
      sharp_rc->commitAcceptedTimestepTransportHistory();
}

void
ReducedPressurePIMPLESolve::advanceOuterIterationHistories()
{
  if (_should_solve_pressure)
    advancePressureOuterIterationHistory();

  if (_should_solve_momentum)
    // Keep the full nonlinear history on the previous outer-corrector state
    // for the whole current outer loop. The stock momentum solve shifts this
    // stack every predictor solve; here we only want to advance it
    // once per outer SIMPLE iteration.
    advanceSystemOuterIterationHistory(_momentum_systems);

  // Keep the true timestep-old alpha in solutionOld(), but advance the
  // nonlinear-state stack once per outer iteration so we have a separate
  // previous-outer iterate available.
  if (solvesVolumeFraction())
    advanceSystemOuterIterationHistory(_active_scalar_systems);
}

void
ReducedPressurePIMPLESolve::solveVolumeFractionBeforeFlowCorrection(
    ResidualStorage & residual_storage, const SolverParams & solver_params)
{
  prepareVOFTransportStateForOuterIteration();

  _problem.execute(EXEC_NONLINEAR);
  Moose::PetscSupport::petscSetOptions(_active_scalar_petsc_options, solver_params);
  const auto vf_residuals = solveVolumeFractionSystems();

  adoptPublishedVOFTransportState();

  _problem.execute(EXEC_NONLINEAR);
  storeActiveScalarResiduals(residual_storage, vf_residuals);
}

void
ReducedPressurePIMPLESolve::prepareVOFTransportStateForOuterIteration() const
{
  const bool use_previous_timestep_transport_flux =
      _problem.timeStep() == 1 && _current_outer_iteration == 1;

  if (auto * sharp_rc = sharpInterfaceRC())
  {
    sharp_rc->clearVOFTransportState();
    sharp_rc->freezeVOFTransportState(use_previous_timestep_transport_flux);
  }
}

void
ReducedPressurePIMPLESolve::adoptPublishedVOFTransportState() const
{
  if (auto * sharp_rc = sharpInterfaceRC())
    sharp_rc->adoptPublishedVOFTransportState();
}

void
ReducedPressurePIMPLESolve::storeActiveScalarResiduals(
    ResidualStorage & residual_storage,
    const std::vector<std::pair<unsigned int, Real>> & vf_residuals) const
{
  for (const auto i : index_range(vf_residuals))
    residual_storage.ns_residuals[residual_storage.active_scalar_indices[i]] = vf_residuals[i];
}

bool
ReducedPressurePIMPLESolve::shouldAssembleMomentumPredictorWithoutSolve() const
{
  return _should_solve_pressure && !_momentum_systems.empty() && _rc_uo;
}

bool
ReducedPressurePIMPLESolve::shouldSolveEnergyAfterFlowLoop() const
{
  return !_has_low_mach_mass_system;
}

bool
ReducedPressurePIMPLESolve::shouldSolveActiveScalarsAfterFlowLoop() const
{
  return false;
}

void
ReducedPressurePIMPLESolve::finalizeSolve(const bool converged)
{
  if (auto * sharp_rc = sharpInterfaceRC())
    sharp_rc->setSuppressExplicitHydrostaticPressureFlux(false);

  if (converged && _low_mach_mass_flux)
  {
    _low_mach_mass_flux->synchronizeDensityToEOS();
    _problem.execute(EXEC_NONLINEAR);
  }
}

void
ReducedPressurePIMPLESolve::addMomentumPredictorExplicitForcing(const unsigned int system_i,
                                                                NumericVector<Number> & rhs)
{
  if (auto * sharp_rc = sharpInterfaceRC(); sharp_rc && sharp_rc->splitMomentumPredictorOperator())
    sharp_rc->addMomentumPredictorExplicitForcing(system_i, rhs);
}

ConservativeSharpInterfaceRhieChowMassFlux *
ReducedPressurePIMPLESolve::sharpInterfaceRC() const
{
  return dynamic_cast<ConservativeSharpInterfaceRhieChowMassFlux *>(_rc_uo);
}

void
ReducedPressurePIMPLESolve::commitAcceptedTimestepTransportHistory() const
{
  if (auto * sharp_rc = sharpInterfaceRC())
    sharp_rc->commitAcceptedTimestepTransportHistory();
}

void
ReducedPressurePIMPLESolve::synchronizeSystemState(LinearSystem & system) const
{
  auto & current_local_solution = *(system.system().current_local_solution);
  current_local_solution.close();
  system.setSolution(current_local_solution);

  auto & current_solution = system.solution();
  current_solution.close();

  system.solutionOld().close();
  system.solutionOld() = current_solution;
  system.solutionOld().close();

  for (unsigned int state = 1;
       system.hasSolutionState(state, Moose::SolutionIterationType::Nonlinear);
       ++state)
  {
    auto & nonlinear_state = system.solutionState(state, Moose::SolutionIterationType::Nonlinear);
    nonlinear_state.close();
    nonlinear_state = system.solutionOld();
    nonlinear_state.close();
  }

  if (auto * previous_newton_solution = system.solutionPreviousNewton())
  {
    previous_newton_solution->close();
    *previous_newton_solution = current_solution;
    previous_newton_solution->close();
  }
}

void
ReducedPressurePIMPLESolve::setPreviousNewtonToCurrent(LinearSystem & system) const
{
  if (auto * previous_solution = system.solutionPreviousNewton())
  {
    *previous_solution = *(system.system().current_local_solution);
    previous_solution->close();
  }
}

void
ReducedPressurePIMPLESolve::advanceVolumeFractionSubcycleOldState(LinearSystem & system) const
{
  system.solutionOld() = *(system.system().current_local_solution);
  system.solutionOld().close();

  if (auto * previous_solution = system.solutionPreviousNewton())
  {
    *previous_solution = system.solutionOld();
    previous_solution->close();
  }
}

void
ReducedPressurePIMPLESolve::initializeStartupPressureField(const SolverParams & solver_params)
{
  if (!startupPressureInitializationEnabled() || _problem.timeStep() != 1)
    return;

  if (!_should_solve_pressure)
    return;

  _console << "Applying startup continuity projection before PIMPLE iterations" << std::endl;

  // Honor the current reduced-pressure field, assemble the momentum predictor coefficients, and run
  // pressure-only startup continuity corrections before the first outer iteration.
  if (_momentum_systems.empty() || !_rc_uo)
    mooseError("ReducedPressurePIMPLESolve startup projection requires momentum systems and a "
               "Rhie-Chow user object.");

  assembleMomentumPredictorWithoutSolve();
  _rc_uo->initFaceMassFlux();

  _console << "Applying startup continuity corrections" << std::endl;

  for (unsigned int startup_it = 0; startup_it < _startup_flux_corrections; ++startup_it)
    (void)correctStartupContinuityOnce(solver_params);

  synchronizeSystemState(_pressure_system);
  _problem.execute(EXEC_NONLINEAR);
}

void
ReducedPressurePIMPLESolve::postPreparePressureCorrectorState(const bool subtract_updated_pressure)
{
  if (auto * sharp_rc = sharpInterfaceRC())
    sharp_rc->updateAdditionalPressureFluxFunctors(subtract_updated_pressure, _print_fields);

  // Refresh the patch velocity / target-flux state from the latest momentum
  // predictor before assembling the constrained pressure boundary gradient.
  _rc_uo->updateVelocityBoundaryState();

  _rc_uo->updatePressureBoundaryNormalGradients(/* apply_pressure_flux_adjustment = */ false);
}

unsigned int
ReducedPressurePIMPLESolve::computeVolumeFractionSubcycles() const
{
  unsigned int subcycles = _volume_fraction_subcycles;

  if (const auto * sharp_rc = sharpInterfaceRC())
  {
    const Real alpha_courant = sharp_rc->maxCourant(_problem.dt());
    if (std::isfinite(alpha_courant) && alpha_courant > _volume_fraction_max_courant)
    {
      const auto required_subcycles =
          static_cast<unsigned int>(std::ceil(alpha_courant / _volume_fraction_max_courant));
      subcycles = std::max(subcycles, std::max(required_subcycles, 1u));
    }
  }

  return std::max(subcycles, 1u);
}

void
ReducedPressurePIMPLESolve::setProblemSubcycleTime(const unsigned int subcycle,
                                                   const Real subcycle_dt,
                                                   const Real global_time_old)
{
  _problem.dt() = subcycle_dt;
  _problem.timeOld() = global_time_old + subcycle * subcycle_dt;
  _problem.time() = _problem.timeOld() + subcycle_dt;
}

std::vector<std::pair<unsigned int, Real>>
ReducedPressurePIMPLESolve::solveVolumeFractionSystems()
{
  ProblemTimeGuard time_guard(_problem);

  std::vector<std::pair<unsigned int, Real>> residuals(_active_scalar_system_names.size(),
                                                       std::make_pair(0, 1.0));

  const Real global_dt = _problem.dt();
  const Real global_time = _problem.time();
  const Real global_time_old = _problem.timeOld();
  const unsigned int num_subcycles = computeVolumeFractionSubcycles();
  const Real subcycle_dt = global_dt / num_subcycles;

  if (num_subcycles > _volume_fraction_subcycles)
    _console << name() << ": increasing alpha subcycles from " << _volume_fraction_subcycles
             << " to " << num_subcycles << " to keep alpha CFL <= " << _volume_fraction_max_courant
             << " at dt=" << global_dt << std::endl;

  for (const auto i : index_range(_active_scalar_system_names))
    residuals[i] =
        solveOneVolumeFractionSystem(i, num_subcycles, subcycle_dt, global_dt, global_time_old);

  _problem.dt() = global_dt;
  _problem.time() = global_time;
  _problem.timeOld() = global_time_old;
  finalizeVolumeFractionTransportState();

  return residuals;
}

std::pair<unsigned int, Real>
ReducedPressurePIMPLESolve::solveOneVolumeFractionSystem(const unsigned int i,
                                                         const unsigned int num_subcycles,
                                                         const Real subcycle_dt,
                                                         const Real global_dt,
                                                         const Real global_time_old)
{
  auto * system = _active_scalar_systems[i];
  auto * corrector = sharpInterfaceVOFCorrector(_active_scalar_system_names[i]);

  if (!corrector)
  {
    auto * level_set_corrector = massConservedLevelSetCorrector(_active_scalar_system_names[i]);
    if (!activeScalarUsesDeferredCorrection(_active_scalar_system_numbers[i]))
      mooseError("ReducedPressurePIMPLESolve requires either a "
                 "ConservativeSharpInterfaceVOFMULESCorrector or a bounded cubic-upwind "
                 "LinearFVMaterialAdvection kernel for active-scalar system '",
                 _active_scalar_system_names[i],
                 "'.");

    setPreviousNewtonToCurrent(*system);
    auto residual = solveActiveScalarSystem(i);
    if (level_set_corrector)
    {
      // Projection magnitude is not an outer residual: the same nonzero correction can reproduce
      // an unchanged corrected level set. The coupled outer-state snapshot measures convergence.
      level_set_corrector->correct();
      if (_has_low_mach_mass_system)
        applyLevelSetAitkenRelaxation(i, *system, *level_set_corrector);
    }
    return residual;
  }

  if (num_subcycles > 1)
    for (const auto & time_integrator : system->getTimeIntegrators())
      if (time_integrator->numStatesRequired() > 1)
        mooseError("Volume-fraction subcycling is not compatible with the configured multistep "
                   "time integrator in active-scalar system '",
                   _active_scalar_system_names[i],
                   "'. Reduce the global time step to satisfy volume_fraction_max_courant without "
                   "subcycling.");

  system->saveOldSolutions();

  // solutionOld() must stay as the true timestep-old alpha for the whole
  // outer loop. solutionPreviousNewton() is only the local/subcycle field-
  // relaxation state, while the previous-outer iterate now lives in the
  // nonlinear solution-state stack advanced at outer-loop entry.
  setPreviousNewtonToCurrent(*system);

  corrector->resetSubcycleFluxes();

  std::pair<unsigned int, Real> residual{0, 1.0};
  for (const auto subcycle : make_range(num_subcycles))
    residual = runOneVolumeFractionSubcycle(
        i, *system, *corrector, subcycle, subcycle_dt, global_dt, global_time_old);

  system->restoreOldSolutions();
  setPreviousNewtonToCurrent(*system);

  return residual;
}

std::pair<unsigned int, Real>
ReducedPressurePIMPLESolve::runOneVolumeFractionSubcycle(
    const unsigned int i,
    LinearSystem & system,
    ConservativeSharpInterfaceVOFMULESCorrector & corrector,
    const unsigned int subcycle,
    const Real subcycle_dt,
    const Real global_dt,
    const Real global_time_old)
{
  setProblemSubcycleTime(subcycle, subcycle_dt, global_time_old);

  if (subcycle > 0)
    advanceVolumeFractionSubcycleOldState(system);

  // Outer corrections recompute the same physical time step. Only later alpha subcycles advance
  // from the preceding subcycle state.
  corrector.cachePreSubcycleAlpha(subcycle == 0);
  _problem.execute(EXEC_NONLINEAR);

  // This path bounds alpha through the limited face fluxes instead of projecting the solved field.
  // Do not apply the inherited scalar lower limiter here.
  const auto residual = solveAdvectedSystem(_active_scalar_system_numbers[i],
                                            system,
                                            _active_scalar_equation_relaxation[i],
                                            _active_scalar_linear_control,
                                            _active_scalar_l_abs_tol,
                                            1.0,
                                            std::numeric_limits<Real>::min());
  system.computeGradients();
  corrector.applyCorrection(subcycle_dt, subcycle_dt / global_dt);

  return residual;
}

void
ReducedPressurePIMPLESolve::finalizeVolumeFractionTransportState()
{
  for (const auto & system : _active_scalar_systems)
    system->computeGradients();

  for (const auto & system_name : _active_scalar_system_names)
    if (auto * corrector = sharpInterfaceVOFCorrector(system_name))
      corrector->refreshPublishedRhoPhi();
}

std::pair<unsigned int, Real>
ReducedPressurePIMPLESolve::correctStartupContinuityOnce(const SolverParams & solver_params)
{
  PressureStateGuard pressure_state(_pressure_system);
  SharpInterfaceStartupProjectionGuard projection_scope(sharpInterfaceRC());

  preparePressureCorrectorState(true);

  const auto residuals = applyPressureCorrectionStage(true, false, solver_params);

  // Restore the user/equilibrium startup reduced-pressure field. Startup
  // continuity cleanup should repair phi, not overwrite the physical p_rgh field
  // before the first real pressure equation.
  pressure_state.restore();
  _pressure_system.computeGradients();

  return residuals;
}

void
ReducedPressurePIMPLESolve::postPublishPressureCorrectedState()
{
  _rc_uo->updateVelocityBoundaryState();
}
