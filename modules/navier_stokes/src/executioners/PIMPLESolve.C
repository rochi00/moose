//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PIMPLESolve.h"
#include "FEProblem.h"
#include "SegregatedSolverUtils.h"
#include "LinearSystem.h"

#include <limits>

using namespace libMesh;

InputParameters
PIMPLESolve::validParams()
{
  InputParameters params = LinearAssemblySegregatedSolve::validParams();
  params.addParam<unsigned int>(
      "num_piso_iterations",
      0,
      "The maximum number of additional PISO pressure-correction stages without recomputing the "
      "momentum matrix.");
  params.addParam<unsigned int>(
      "num_pressure_nonorthogonal_correctors",
      0,
      "Number of additional non-final pressure equation solves inside each pressure-corrector "
      "stage. A value of 0 means one final pressure solve, 1 means one non-final solve followed "
      "by one final solve, etc.");

  return params;
}

PIMPLESolve::PIMPLESolve(Executioner & ex)
  : LinearAssemblySegregatedSolve(ex),
    _num_piso_iterations(getParam<unsigned int>("num_piso_iterations")),
    _num_pressure_nonorthogonal_correctors(
        getParam<unsigned int>("num_pressure_nonorthogonal_correctors"))
{
}

std::pair<unsigned int, Real>
PIMPLESolve::correctVelocity(const bool /*recompute_face_mass_flux*/,
                             const SolverParams & solver_params)
{
  std::pair<unsigned int, Real> residual;
  unsigned int piso_iteration_counter = 0;
  const bool reconstructed = _rc_uo && _rc_uo->usingReconstructedPressureGradientMethod();
  while (piso_iteration_counter <= _num_piso_iterations)
  {
    const bool last_piso_corrector = piso_iteration_counter == _num_piso_iterations;

    // Reconstructed candidates always need the face flux produced by their own pressure
    // solve, so every corrector must recompute it. Ordinary (non-reconstructed) PIMPLE only
    // needs the final corrector's flux, which feeds the advection terms for the next outer
    // iteration/time step.
    const bool recompute_flux = reconstructed || last_piso_corrector;

    residual = LinearAssemblySegregatedSolve::correctVelocity(recompute_flux, solver_params);

    // After each PISO corrector except the last, refresh only the lagged velocity gradient. H/A
    // continues using the coupling-gradient snapshot associated with the existing momentum system.
    if (!last_piso_corrector && reconstructed)
      _rc_uo->preparePISOCorrector();

    piso_iteration_counter++;
  }

  _current_piso_iteration = 0;

  return residual;
}

void
PIMPLESolve::postPressureCorrectorSolve(const bool /*final_nonorthogonal_corrector*/)
{
  // solvePressureCorrector has already installed the solution and cached the pressure-equation
  // flux. What the next non-orthogonal pass additionally needs is the cell gradients of the
  // pressure it just produced, since those carry the explicit correction.
  _pressure_system.computeGradients();
}

void
PIMPLESolve::advanceSystemOuterIterationHistory(const std::vector<LinearSystem *> & systems) const
{
  for (auto * system : systems)
  {
    unsigned int max_state = 0;
    while (system->hasSolutionState(max_state + 1, Moose::SolutionIterationType::Nonlinear))
      ++max_state;

    for (unsigned int state = max_state; state > 1; --state)
    {
      auto & nonlinear_state =
          system->solutionState(state, Moose::SolutionIterationType::Nonlinear);
      nonlinear_state = system->solutionState(state - 1, Moose::SolutionIterationType::Nonlinear);
      nonlinear_state.close();
    }

    if (max_state >= 1)
    {
      auto & previous_outer_solution =
          system->solutionState(1, Moose::SolutionIterationType::Nonlinear);
      previous_outer_solution = *(system->system().current_local_solution);
      previous_outer_solution.close();
    }
  }
}

void
PIMPLESolve::advancePressureOuterIterationHistory() const
{
  advanceSystemOuterIterationHistory({&_pressure_system});
}
