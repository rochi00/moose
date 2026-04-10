//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVFluxKernel.h"
#include "RhieChowMassFluxMultiPhase.h"
#include "LinearFVAdvectionDiffusionBC.h"

/**
 * Advection kernel for the level-set variable in a CLSVOF solver.
 * Uses the volumetric face flux from RhieChowMassFluxMultiPhase
 * with standard upwind/TVD interpolation. No MULES, no compression.
 */
class LinearFVLevelSetAdvection : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();
  LinearFVLevelSetAdvection(const InputParameters & params);

  virtual Real computeElemMatrixContribution() override;
  virtual Real computeNeighborMatrixContribution() override;
  virtual Real computeElemRightHandSideContribution() override;
  virtual Real computeNeighborRightHandSideContribution() override;
  virtual Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;
  virtual Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;
  virtual void setupFaceData(const FaceInfo * face_info) override;

protected:
  /// The Rhie-Chow user object that provides the volumetric face flux
  const RhieChowMassFluxMultiPhase & _mass_flux_provider;

private:
  /// Cached interpolation coefficients on the current face
  std::pair<Real, Real> _advected_interp_coeffs;

  /// Cached volumetric face flux
  Real _volumetric_face_flux;

  /// The interpolation method to use for the advected quantity
  Moose::FV::InterpMethod _advected_interp_method;
};
