//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVElementalKernel.h"
#include "FaceArgInterface.h"

/**
 * Kernel that adds the surface tension force component to the momentum equation RHS.
 * Computes F = -sigma * kappa * grad(alpha), where kappa is the interface curvature
 * and grad(alpha) localizes the force to the interface.
 *
 * When an external curvature variable is provided (e.g. from a level-set field),
 * that precomputed curvature is used instead of computing it inline from alpha gradients.
 */
class LinearFVMomentumSurfaceTensionForce : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  LinearFVMomentumSurfaceTensionForce(const InputParameters & params);

  virtual Real computeMatrixContribution() override;

  virtual Real computeRightHandSideContribution() override;

protected:

  MooseLinearVariableFV<Real> & getAlphaVariable(const std::string & vname);

  /// The dimension of the problem
  const unsigned int _dim;

  /// Index x|y|z of the momentum equation component
  const unsigned int _index;

  /// The surface tension value
  const Moose::Functor<Real> & _sigma;

  /// The phase fraction (used for grad(alpha) localization and inline curvature)
  MooseLinearVariableFV<Real> & _alpha;

  /// Optional external curvature functor (e.g. from level-set). When provided,
  /// skips inline curvature computation from alpha.
  const Moose::Functor<Real> * const _curvature_functor;
};
