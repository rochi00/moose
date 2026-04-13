//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PLICReconstruction.h"
#include "MooseLinearVariableFV.h"
#include "FEProblemBase.h"
#include "LinearSystem.h"
#include "MooseMesh.h"
#include "ElemInfo.h"
#include "PLICUtils.h"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/elem.h"

registerMooseObject("NavierStokesApp", PLICReconstruction);

InputParameters
PLICReconstruction::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Computes PLIC interface reconstruction for each interfacial cell. "
      "The plane normal comes from grad(phi) and the offset is determined "
      "by matching the cell volume fraction alpha. 2D QUAD4 only.");
  params.addRequiredParam<VariableName>("alpha_variable", "The VOF phase fraction variable.");
  params.addParam<VariableName>("level_set_variable",
                                "The level-set variable (gradient provides normal). "
                                "If not provided, grad(alpha) is used instead (Youngs method).");
  params.addParam<Real>(
      "alpha_tolerance", 1e-6, "Cells with alpha in (tol, 1-tol) are reconstructed.");
  return params;
}

PLICReconstruction::PLICReconstruction(const InputParameters & params)
  : GeneralUserObject(params),
    _alpha_var(dynamic_cast<MooseLinearVariableFV<Real> &>(
        _fe_problem.getVariable(0, getParam<VariableName>("alpha_variable")))),
    _alpha_system(_fe_problem.getLinearSystem(_alpha_var.sys().number())),
    _use_level_set(isParamValid("level_set_variable")),
    _phi_var(_use_level_set ? &dynamic_cast<MooseLinearVariableFV<Real> &>(
                                  _fe_problem.getVariable(
                                      0, getParam<VariableName>("level_set_variable")))
                            : nullptr),
    _phi_system(_use_level_set ? &_fe_problem.getLinearSystem(_phi_var->sys().number()) : nullptr),
    _alpha_tol(getParam<Real>("alpha_tolerance"))
{
  if (_use_level_set)
    _phi_var->computeCellGradients();
  else
    _alpha_var.computeCellGradients();
}

bool
PLICReconstruction::hasPlane(dof_id_type elem_id) const
{
  return _planes.count(elem_id) > 0;
}

const PLICReconstruction::PLICPlane &
PLICReconstruction::getPlane(dof_id_type elem_id) const
{
  return _planes.at(elem_id);
}

void
PLICReconstruction::execute()
{
  auto & mesh = _fe_problem.mesh();
  _planes.clear();

  // Ensure gradients are up to date
  if (_use_level_set)
    _phi_system->computeGradients();
  else
    _alpha_system.computeGradients();

  // Get alpha solution vector
  const auto & alpha_sys =
      libMesh::cast_ref<libMesh::LinearImplicitSystem &>(_alpha_var.sys().system());
  const auto & alpha_solution = *alpha_sys.solution;
  const auto alpha_sys_num = _alpha_var.sys().number();
  const auto alpha_var_num = _alpha_var.number();

  unsigned int n_interface = 0, n_zero_grad = 0, n_total = 0;
  Real alpha_min = 1.0, alpha_max = 0.0;

  for (const auto * elem_info_ptr : mesh.elemInfoVector())
  {
    const auto & ei = *elem_info_ptr;
    const auto * elem = ei.elem();
    n_total++;

    // Read alpha for this cell
    const auto alpha_dof = ei.dofIndices()[alpha_sys_num][alpha_var_num];
    const Real alpha = alpha_solution(alpha_dof);

    if (alpha < alpha_min) alpha_min = alpha;
    if (alpha > alpha_max) alpha_max = alpha;

    // Skip pure cells
    if (alpha <= _alpha_tol || alpha >= 1.0 - _alpha_tol)
      continue;

    n_interface++;

    // Get interface normal from grad(phi) or grad(alpha)
    const auto grad_field = _use_level_set ? _phi_var->gradSln(ei) : _alpha_var.gradSln(ei);
    const Real grad_mag = grad_field.norm();

    if (grad_mag < 1.0e-14)
    {
      n_zero_grad++;
      continue;
    }

    VectorValue<Real> n_hat = grad_field / grad_mag;

    // Get cell dimensions (assumes axis-aligned QUAD4)
    const Point & p0 = elem->point(0);
    const Point & p2 = elem->point(2); // diagonally opposite corner
    const Real dx = std::abs(p2(0) - p0(0));
    const Real dy = std::abs(p2(1) - p0(1));

    // Lower-left corner
    const Point x0(std::min(p0(0), p2(0)), std::min(p0(1), p2(1)), 0.0);

    // Compute local PLIC offset and convert to global
    const Real d_local = NS::PLIC::computePlaneOffset2DRect(n_hat, alpha, dx, dy);
    const Real d_global = NS::PLIC::localToGlobalOffset(n_hat, d_local, x0, dx, dy);

    _planes[elem->id()] = {n_hat, d_global};
  }

  mooseWarning("PLICReconstruction: ",
              n_total, " total, ",
              n_interface, " interface, ",
              _planes.size(), " planes, ",
              "alpha=[", alpha_min, ",", alpha_max, "]");
}
