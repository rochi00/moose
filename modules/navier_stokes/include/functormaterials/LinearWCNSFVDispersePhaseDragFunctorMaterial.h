//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"
#include "NonADFunctorInterface.h"

/**
 * Computes the drag coefficient between the dispersed and the continuous phase of a two-phase
 * mixture, for the linear finite volume discretization.
 *
 * This is the linear finite volume counterpart of NSFVDispersePhaseDragFunctorMaterial. It computes
 * in plain reals, since the linear finite volume discretization assembles a matrix directly and has
 * no use for the derivatives the nonlinear version carries.
 */
class LinearWCNSFVDispersePhaseDragFunctorMaterial : public FunctorMaterial,
                                                     public NonADFunctorInterface
{
public:
  static InputParameters validParams();
  LinearWCNSFVDispersePhaseDragFunctorMaterial(const InputParameters & parameters);

protected:
  /// The dimension of the simulation
  const unsigned int _dim;

  /// Slip velocity in the x direction
  const Moose::Functor<Real> & _u_var;
  /// Velocity in the y direction
  const Moose::Functor<Real> * const _v_var;
  /// Velocity in the z direction
  const Moose::Functor<Real> * const _w_var;

  /// Density used to form the particle Reynolds number
  const Moose::Functor<Real> & _rho;

  /// Dynamic viscosity used to form the particle Reynolds number
  const Moose::Functor<Real> & _mu;

  /// Particle diameter in the dispersed phase
  const Moose::Functor<Real> & _particle_diameter;
};
