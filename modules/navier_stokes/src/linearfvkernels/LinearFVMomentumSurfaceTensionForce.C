//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVMomentumSurfaceTensionForce.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "NS.h"
#include "FEProblemBase.h"

registerMooseObject("NavierStokesApp", LinearFVMomentumSurfaceTensionForce);

InputParameters
LinearFVMomentumSurfaceTensionForce::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription("Volumetric force imposed by the surface tension.");
  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>(
      "momentum_component",
      momentum_component,
      "The component of the momentum equation that this kernel applies to.");
  params.addRequiredParam<MooseFunctorName>("sigma", "The value of the surface tension.");
  params.addRequiredParam<VariableName>("alpha", "The phase fraction");
  params.addParam<MooseFunctorName>(
      "curvature",
      "Optional precomputed curvature (e.g. from a level-set field). When provided, "
      "this value is used instead of computing curvature inline from alpha gradients.");
  return params;
}

LinearFVMomentumSurfaceTensionForce::LinearFVMomentumSurfaceTensionForce(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _index(getParam<MooseEnum>("momentum_component")),
    _sigma(getFunctor<Real>("sigma")),
    _alpha(getAlphaVariable("alpha")),
    _curvature_functor(isParamValid("curvature") ? &getFunctor<Real>("curvature") : nullptr)
{
  _alpha.computeCellGradients();
}

MooseLinearVariableFV<Real> &
LinearFVMomentumSurfaceTensionForce::getAlphaVariable(const std::string & vname)
{
  auto * ptr = dynamic_cast<MooseLinearVariableFV<Real> *>(
      &_fe_problem.getVariable(_tid, getParam<VariableName>(vname)));

  if (!ptr)
    paramError("alpha", "The alpha variable should be of type MooseLinearVariableFVReal!");

  return *ptr;
}

Real
LinearFVMomentumSurfaceTensionForce::computeMatrixContribution()
{
  return 0.0;
}

Real
LinearFVMomentumSurfaceTensionForce::computeRightHandSideContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const auto coord_type = _subproblem.mesh().getCoordSystem(elem_arg.elem->subdomain_id());
  const auto rz_radial_coord = _subproblem.mesh().getAxisymmetricRadialCoord();

  Real curvature = 0.0;

  std::vector<RealVectorValue> A_sys;
  std::vector<Real> rhs_sys;
  unsigned int nrow = 0;

  auto action_functor = [this,
                         &elem_arg,
                         &curvature,
                         &A_sys,
                         &rhs_sys,
                         &nrow](const Elem & elem,
                                const Elem * /*neighbor*/,
                                const FaceInfo * const fi,
                                const Point & surface_vector,
                                Real /*coord*/,
                                const bool /*elem_has_info*/)
  {
      mooseAssert(fi, "We need a FaceInfo for this action_functor");

      Moose::FaceArg face_arg = Moose::FaceArg{fi,
                                               Moose::FV::LimiterType::CentralDifference,
                                               true,
                                               /* correct_skewness */ false,
                                               this->_current_elem_info->elem(),
                                               nullptr};

      // Compute inline curvature from alpha gradients only when no external curvature is provided
      if (!this->_curvature_functor)
      {
        const auto grad_alpha = MetaPhysicL::raw_value(this->_alpha.gradient(face_arg, this->determineState()));
        constexpr Real tiny = 1.0e-14;
        const auto grad_alpha_norm = grad_alpha / (grad_alpha.norm() + tiny);
        curvature += grad_alpha_norm * surface_vector;
      }

      // Always reconstruct grad(alpha) via SVD for the force localization
      A_sys.push_back(fi->faceCentroid() - elem.vertex_average());

      const auto alpha_face = MetaPhysicL::raw_value(_alpha(face_arg, determineState()));
      const auto alpha_elem = MetaPhysicL::raw_value(_alpha(elem_arg, determineState()));
      rhs_sys.push_back(alpha_face - alpha_elem);

      ++nrow;
  };

  Moose::FV::loopOverElemFaceInfo(
      *_current_elem_info->elem(), _subproblem.mesh(), action_functor, coord_type, rz_radial_coord);

  // Use external curvature if provided, otherwise use inline computation
  if (_curvature_functor)
    curvature = (*_curvature_functor)(elem_arg, determineState());
  else
    curvature = curvature / _current_elem_volume;

  // Reconstruct grad(alpha) at cell center via least-squares (SVD)
  DenseMatrix<Real> A(nrow, _dim);
  DenseVector<Real> b(nrow), x(_dim);

  for (unsigned int i = 0; i < nrow; ++i)
  {
    for (unsigned int j = 0; j < _dim; ++j)
      A(i, j) = A_sys[i](j);
    b(i) = rhs_sys[i];
  }

  A.svd_solve(b, x);

  return -_sigma(elem_arg, determineState()) * curvature * x(_index) * _current_elem_volume;
}