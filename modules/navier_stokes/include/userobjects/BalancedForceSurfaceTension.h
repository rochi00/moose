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
#include "NonADFunctorInterface.h"

/**
 * Provides face-interpolated curvature and cell alpha values for the
 * balanced-force Rhie-Chow surface tension correction.
 *
 * The balanced-force approach adds a surface tension flux to the face
 * mass flux using the SAME discrete gradient operator as the pressure:
 *
 *   st_flux = sigma * kappa_f * (alpha_N * mc_N + alpha_C * mc_C)
 *
 * where mc_N, mc_C are the pressure diffusion kernel's matrix contributions.
 * Since (alpha_N - alpha_C) = 0 far from the interface, garbage curvature
 * values in pure cells are harmlessly multiplied by zero.
 *
 * At equilibrium: grad(p) = sigma * kappa * grad(alpha), so
 * p_grad_flux = st_flux and the face mass flux is zero.
 */
class BalancedForceSurfaceTension : public GeneralUserObject, public NonADFunctorInterface
{
public:
  static InputParameters validParams();
  BalancedForceSurfaceTension(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual void finalize() override {}

  /// Get the curvature at a face (interpolated from cell-centered values)
  Real getFaceCurvature(const FaceInfo & fi) const;

  /// Get the alpha (VOF fraction) value for a cell
  Real getAlpha(dof_id_type elem_id) const;

  /// Get the surface tension coefficient
  Real sigma() const { return _sigma; }

protected:
  /// Surface tension coefficient
  const Real _sigma;

  /// Curvature functor (reads from AuxVariable populated by LinearFVCurvatureAux)
  const Moose::Functor<Real> & _kappa;

  /// VOF phase fraction variable
  MooseLinearVariableFV<Real> & _alpha_var;
};
