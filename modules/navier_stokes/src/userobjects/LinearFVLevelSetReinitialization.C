//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVLevelSetReinitialization.h"
#include "MooseLinearVariableFV.h"
#include "FEProblemBase.h"
#include "LinearSystem.h"
#include "MooseMesh.h"
#include "ElemInfo.h"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/petsc_vector.h"

registerMooseObject("NavierStokesApp", LinearFVLevelSetReinitialization);

InputParameters
LinearFVLevelSetReinitialization::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Performs Sussman PDE-based reinitialization on a level-set variable to "
      "restore the signed distance property |grad(phi)| = 1.");
  params.addRequiredParam<VariableName>("level_set_variable",
                                        "The level-set variable to reinitialize.");
  params.addParam<unsigned int>(
      "num_iterations", 5, "Number of pseudo-time iterations for reinitialization.");
  params.addParam<Real>("cfl", 0.5, "CFL number for the pseudo-time step (dt_pseudo = cfl * h).");
  params.addParam<Real>(
      "epsilon_factor",
      1.5,
      "The smoothing parameter epsilon = epsilon_factor * h_min, where h_min is the "
      "minimum cell size. Controls the width of the sign function transition.");
  return params;
}

LinearFVLevelSetReinitialization::LinearFVLevelSetReinitialization(const InputParameters & params)
  : GeneralUserObject(params),
    _phi_var(dynamic_cast<MooseLinearVariableFV<Real> &>(
        _fe_problem.getVariable(0, getParam<VariableName>("level_set_variable")))),
    _phi_system(_fe_problem.getLinearSystem(_phi_var.sys().number())),
    _phi_implicit_system(
        libMesh::cast_ref<libMesh::LinearImplicitSystem &>(_phi_system.system())),
    _num_iterations(getParam<unsigned int>("num_iterations")),
    _cfl(getParam<Real>("cfl")),
    _epsilon_factor(getParam<Real>("epsilon_factor")),
    _dt_pseudo(0.0),
    _epsilon(0.0)
{
  _phi_var.computeCellGradients();
}

void
LinearFVLevelSetReinitialization::execute()
{
  if (_num_iterations == 0)
    return;

  auto & mesh = _fe_problem.mesh();

  // Compute h_min, dt_pseudo, and epsilon (mesh doesn't change, but do it once per execute
  // in case of AMR or mesh changes between timesteps)
  Real h_min = std::numeric_limits<Real>::max();
  for (const auto * elem_info_ptr : mesh.elemInfoVector())
    if (elem_info_ptr->volume() > 0)
    {
      const Real h = std::pow(elem_info_ptr->volume(), 1.0 / mesh.dimension());
      h_min = std::min(h_min, h);
    }
  _communicator.min(h_min);
  _dt_pseudo = _cfl * h_min;
  _epsilon = _epsilon_factor * h_min;

  // Get the current solution vector
  NumericVector<Number> & solution = *_phi_implicit_system.solution;

  // Store the original φ₀ (before reinitialization) for the sign function
  auto phi_0 = solution.clone();
  auto & petsc_phi_0 = dynamic_cast<PetscVector<Number> &>(*phi_0);

  // Allocate reusable buffer for the explicit update
  if (!_phi_old)
    _phi_old = solution.clone();

  auto & petsc_solution = dynamic_cast<PetscVector<Number> &>(solution);
  auto & petsc_phi_old = dynamic_cast<PetscVector<Number> &>(*_phi_old);

  const auto sys_num = _phi_var.sys().number();
  const auto var_num = _phi_var.number();

  // Get phi_0 array once outside the iteration loop
  const PetscScalar * phi0_array = petsc_phi_0.get_array_read();

  for (unsigned int iter = 0; iter < _num_iterations; ++iter)
  {
    // Recompute gradients from the current solution
    _phi_system.computeGradients();

    // Copy current solution into reusable buffer (read old, write new)
    *_phi_old = solution;

    PetscScalar * sol_array = petsc_solution.get_array();
    const PetscScalar * old_array = petsc_phi_old.get_array_read();

    for (const auto * elem_info_ptr : mesh.elemInfoVector())
    {
      const auto & ei = *elem_info_ptr;
      const auto dof_id = ei.dofIndices()[sys_num][var_num];

      // Only update locally owned DOFs
      if (dof_id < petsc_solution.first_local_index() ||
          dof_id >= petsc_solution.last_local_index())
        continue;

      const auto local_id = dof_id - petsc_solution.first_local_index();

      // Smoothed sign function: S_ε(φ₀) = φ₀ / sqrt(φ₀² + ε²)
      const Real phi0_val = phi0_array[local_id];
      const Real sign_phi = phi0_val / std::sqrt(phi0_val * phi0_val + _epsilon * _epsilon);

      // Get gradient magnitude |∇φ| at cell center
      const auto grad_phi = _phi_var.gradSln(ei);
      const Real grad_mag = grad_phi.norm();

      // Explicit Euler update: φ^{k+1} = φ^k - Δτ * S(φ₀) * (|∇φ| - 1)
      sol_array[local_id] = old_array[local_id] - _dt_pseudo * sign_phi * (grad_mag - 1.0);
    }

    petsc_phi_old.restore_array();
    petsc_solution.restore_array();

    solution.close();
  }

  // Restore phi_0 read array
  petsc_phi_0.restore_array();

  // Final gradient recomputation for downstream consumers
  _phi_system.computeGradients();
}
