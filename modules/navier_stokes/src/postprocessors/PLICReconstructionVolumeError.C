//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PLICReconstructionVolumeError.h"

#include "PLICReconstruction.h"
#include "PLICUtils.h"
#include "PLICUtils3D.h"
#include "ElemInfo.h"
#include "FEProblemBase.h"
#include "LinearSystem.h"
#include "MooseMesh.h"

#include "libmesh/elem.h"
#include "libmesh/linear_implicit_system.h"

registerMooseObject("NavierStokesApp", PLICReconstructionVolumeError);

InputParameters
PLICReconstructionVolumeError::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Computes the maximum alpha reconstruction mismatch from the planes stored in a "
      "PLICReconstruction user object.");
  params.addRequiredParam<UserObjectName>("plic_uo", "The PLIC reconstruction user object.");
  params.addRequiredParam<VariableName>("alpha_variable", "The phase fraction variable.");
  return params;
}

PLICReconstructionVolumeError::PLICReconstructionVolumeError(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _plic(getUserObject<PLICReconstruction>("plic_uo")),
    _alpha_var(dynamic_cast<MooseLinearVariableFV<Real> &>(
        _fe_problem.getVariable(0, getParam<VariableName>("alpha_variable")))),
    _alpha_system(_fe_problem.getLinearSystem(_alpha_var.sys().number())),
    _max_alpha_error(0.0)
{
}

void
PLICReconstructionVolumeError::initialize()
{
  _max_alpha_error = 0.0;

  const auto & alpha_sys =
      libMesh::cast_ref<libMesh::LinearImplicitSystem &>(_alpha_var.sys().system());
  const auto & alpha_solution = *alpha_sys.solution;
  const auto alpha_sys_num = _alpha_var.sys().number();
  const auto alpha_var_num = _alpha_var.number();

  dof_id_type local_plane_count = 0;

  for (const auto * elem_info_ptr : _fe_problem.mesh().elemInfoVector())
  {
    const auto & ei = *elem_info_ptr;
    const auto * elem = ei.elem();

    if (!_plic.hasPlane(elem->id()))
      continue;

    local_plane_count++;

    const auto alpha_dof = ei.dofIndices()[alpha_sys_num][alpha_var_num];
    const Real alpha = alpha_solution(alpha_dof);
    const auto & plane = _plic.getPlane(elem->id());

    Real reconstructed_alpha = 0.0;
    if (elem->dim() == 2 && elem->type() == QUAD4)
    {
      const Point & p0 = elem->point(0);
      const Point & p2 = elem->point(2);
      const Real dx = std::abs(p2(0) - p0(0));
      const Real dy = std::abs(p2(1) - p0(1));
      const Point x0(std::min(p0(0), p2(0)), std::min(p0(1), p2(1)), 0.0);

      reconstructed_alpha = NS::PLIC::volumeFractionInRect(plane.n_hat, plane.d, x0, dx, dy);
    }
    else
    {
      std::vector<Point> vertices;
      std::vector<std::vector<unsigned int>> faces;
      NS::PLIC::elemGeometry(*elem, vertices, faces);
      reconstructed_alpha =
          NS::PLIC::truncatedVolumeFraction3D(vertices, faces, elem->volume(), plane.n_hat, plane.d);
    }

    _max_alpha_error = std::max(_max_alpha_error, std::abs(reconstructed_alpha - alpha));
  }

  gatherMax(_max_alpha_error);
  gatherSum(local_plane_count);

  if (!local_plane_count)
    mooseError("PLICReconstructionVolumeError expected at least one reconstructed interface plane.");
}

PostprocessorValue
PLICReconstructionVolumeError::getValue() const
{
  return _max_alpha_error;
}
