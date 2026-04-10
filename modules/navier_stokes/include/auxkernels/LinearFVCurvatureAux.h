//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "MooseLinearVariableFV.h"

/**
 * Computes the cell-centred interface curvature κ = ∇·n̂ via the divergence
 * theorem over cell faces, where n̂ = ∇α / |∇α|.  The result is stored in
 * an elemental MooseLinearVariableFVReal aux variable so that other objects
 * (e.g. RhieChowMassFluxMultiPhase, LinearFVMomentumSurfaceTensionForce) can
 * read a consistent κ field through the functor interface.
 */
class LinearFVCurvatureAux : public AuxKernel
{
public:
  static InputParameters validParams();
  LinearFVCurvatureAux(const InputParameters & params);

protected:
  virtual Real computeValue() override;

  MooseLinearVariableFV<Real> & getAlphaVariable(const std::string & vname);

  /// Phase fraction variable whose gradient drives the curvature estimate
  MooseLinearVariableFV<Real> & _alpha;
};
