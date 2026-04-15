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

class PLICReconstruction;

/**
 * Geometric VOF advection kernel using PLIC interface reconstruction.
 *
 * Computes the volume flux of the liquid phase through each face by
 * clipping the advected PLIC plane against the face polygon and
 * integrating in time with Simpson's 3-point quadrature.
 *
 * This is a drop-in replacement for LinearFVMultiPhaseFractionAdvection
 * that eliminates numerical diffusion at the interface.
 *
 * The advection is purely explicit: all contributions go to the RHS.
 * The matrix diagonal comes from LinearFVTimeDerivative (V/dt).
 */
class LinearFVPLICAdvection : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();
  LinearFVPLICAdvection(const InputParameters & params);

  virtual Real computeElemMatrixContribution() override;
  virtual Real computeNeighborMatrixContribution() override;
  virtual Real computeElemRightHandSideContribution() override;
  virtual Real computeNeighborRightHandSideContribution() override;
  virtual Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;
  virtual Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;
  virtual void setupFaceData(const FaceInfo * face_info) override;

private:
  /// The PLIC reconstruction user object (provides interface planes per cell)
  const PLICReconstruction & _plic;

  /// The Rhie-Chow user object (provides volumetric face fluxes)
  const RhieChowMassFluxMultiPhase & _mass_flux_provider;

  /// Cells with alpha in (tol, 1-tol) are treated as interfacial
  const Real _alpha_tol;

  /// Cached geometric flux for the current face (m^3).
  /// Positive means liquid flows from elem to neighbor.
  Real _geometric_flux;

  /// Extract ordered face polygon vertices from a FaceInfo.
  void getFaceVertices(const FaceInfo & fi, std::vector<Point> & verts) const;
};
