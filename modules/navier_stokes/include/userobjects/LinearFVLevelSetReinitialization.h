//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include "MooseLinearVariableFV.h"

class LinearSystem;

namespace libMesh
{
class LinearImplicitSystem;
}

/**
 * Performs PDE-based Sussman reinitialization on a level-set variable to restore
 * the signed distance property |∇φ| = 1.
 *
 * The equation solved in pseudo-time is:
 *   ∂φ/∂τ + S_ε(φ₀)(|∇φ| - 1) = 0
 *
 * where S_ε(φ) = φ / sqrt(φ² + ε²) is a smoothed sign function and φ₀ is the
 * level-set field before reinitialization (preserved to keep the zero contour fixed).
 *
 * Reference: Sussman et al. (1994), "A Level Set Approach for Computing Solutions
 * to Incompressible Two-Phase Flow", J. Comput. Phys. 114, 146-159.
 */
class LinearFVLevelSetReinitialization : public GeneralUserObject
{
public:
  static InputParameters validParams();
  LinearFVLevelSetReinitialization(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}

protected:
  /// The level-set variable to reinitialize
  MooseLinearVariableFV<Real> & _phi_var;

  /// The linear system containing φ
  LinearSystem & _phi_system;

  /// The libMesh implicit system
  libMesh::LinearImplicitSystem & _phi_implicit_system;

  /// Number of pseudo-time iterations
  const unsigned int _num_iterations;

  /// CFL number for the pseudo-time step (Δτ = cfl * h_min)
  const Real _cfl;

  /// Smoothing parameter ε for the sign function (default: h_min)
  const Real _epsilon_factor;

  /// Cached pseudo-time step and epsilon (computed once in initialize)
  Real _dt_pseudo;
  Real _epsilon;

  /// Reusable buffer for the explicit update (avoids cloning each iteration)
  std::unique_ptr<NumericVector<Number>> _phi_old;
};
