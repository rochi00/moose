//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVCurvatureAux.h"
#include "SubProblem.h"
#include "FEProblemBase.h"
#include "MooseMesh.h"
#include "ElemInfo.h"
#include "FVUtils.h"

registerMooseObject("NavierStokesApp", LinearFVCurvatureAux);

InputParameters
LinearFVCurvatureAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Computes cell-centred interface curvature κ = ∇·n̂ via the divergence theorem over cell "
      "faces.  Intended as a pre-processing step for the balanced-force CSF algorithm.");
  params.addRequiredParam<VariableName>("alpha", "The phase fraction variable.");
  return params;
}

LinearFVCurvatureAux::LinearFVCurvatureAux(const InputParameters & params)
  : AuxKernel(params), _alpha(getAlphaVariable("alpha"))
{
  if (isNodal())
    mooseError("LinearFVCurvatureAux only supports elemental (cell-centred) variables.");
  _alpha.computeCellGradients();
}

MooseLinearVariableFV<Real> &
LinearFVCurvatureAux::getAlphaVariable(const std::string & vname)
{
  auto & fep = dynamic_cast<FEProblemBase &>(_subproblem);
  auto * ptr = dynamic_cast<MooseLinearVariableFV<Real> *>(
      &fep.getVariable(_tid, getParam<VariableName>(vname)));
  if (!ptr)
    paramError(vname, "The alpha variable must be of type MooseLinearVariableFVReal!");
  return *ptr;
}

Real
LinearFVCurvatureAux::computeValue()
{
  const auto coord_type = _subproblem.mesh().getCoordSystem(_current_elem->subdomain_id());
  const auto rz_radial_coord = _subproblem.mesh().getAxisymmetricRadialCoord();
  const auto & elem_info = _subproblem.mesh().elemInfo(_current_elem->id());
  const Real volume = elem_info.volume() * elem_info.coordFactor();
  const auto state = determineState();

  Real curvature = 0.0;

  auto action_functor = [this, &curvature, &state](const Elem & /*elem*/,
                                                    const Elem * /*neighbor*/,
                                                    const FaceInfo * const fi,
                                                    const Point & surface_vector,
                                                    Real /*coord*/,
                                                    const bool /*elem_has_info*/)
  {
    mooseAssert(fi, "We need a FaceInfo for this action_functor");
    const Moose::FaceArg face_arg{
        fi, Moose::FV::LimiterType::CentralDifference, true, false, this->_current_elem, nullptr};
    const auto grad_alpha = MetaPhysicL::raw_value(_alpha.gradient(face_arg, state));
    constexpr Real tiny = 1.0e-14;
    curvature += (grad_alpha / (grad_alpha.norm() + tiny)) * surface_vector;
  };

  Moose::FV::loopOverElemFaceInfo(
      *_current_elem, _subproblem.mesh(), action_functor, coord_type, rz_radial_coord);

  return curvature / volume;
}
