//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVInletOutletScalarBC.h"
#include "NSFVUtils.h"
#include "NS.h"

class ElemInfo;

/**
 * Velocity boundary condition of a pressure-controlled inlet/outlet: a zero-gradient outlet on
 * outflow, a fixed tangential velocity with an extrapolated normal velocity on backflow.
 */
class LinearFVPressureInletOutletVelocityBC : public LinearFVInletOutletScalarBC
{
public:
  static InputParameters validParams();

  LinearFVPressureInletOutletVelocityBC(const InputParameters & parameters);

protected:
  Real computeBackflowBoundaryValue() const override;
  Real computeBackflowBoundaryValueMatrixContribution() const override;
  RealVectorValue outwardUnitNormal() const;

  const unsigned int _dim;
  const NS::LinearFVVelocityVariableArray _velocity_vars;
  const unsigned int _index;
};
