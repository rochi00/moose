//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "LinearFVFluxKernel.h"

#include <unordered_map>

class RhieChowMassFlux;

/**
 * Implicit conservative advection operator for the temporary low-Mach density.
 *
 * The matrix uses first-order donor-cell upwinding so the transported density remains bounded.
 * An optional limited material-compression flux may be added explicitly.
 *
 * The exact flux from the last assembly is exposed through faceMassFluxDensity() so the converged
 * discrete mass flux can be published without reconstructing it through a different numerical
 * path.
 */
class LinearFVLowMachMassAdvection : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();

  LinearFVLowMachMassAdvection(const InputParameters & params);

  Real computeElemMatrixContribution() override;
  Real computeNeighborMatrixContribution() override;
  Real computeElemRightHandSideContribution() override;
  Real computeNeighborRightHandSideContribution() override;
  Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;
  Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;
  void setupFaceData(const FaceInfo * face_info) override;

  /**
   * Return the converged mass-flux density, oriented with the FaceInfo normal.
   */
  Real faceMassFluxDensity(const FaceInfo & face_info);

  const MooseFunctorName & volumetricFaceFluxName() const { return _volumetric_face_flux_name; }

protected:
  Real volumetricFaceFlux(const FaceInfo & face_info) const;

  const MooseFunctorName _volumetric_face_flux_name;

  /// Optional Rhie-Chow provider of volumetric face flux.
  const RhieChowMassFlux * const _volumetric_flux_provider;

  /// Optional volumetric face-flux functor, oriented with the FaceInfo normal.
  const Moose::Functor<Real> * const _volumetric_face_flux;

  /// Optional limited artificial material-compression flux and its phase densities.
  const Moose::Functor<Real> * const _material_compression_flux;
  const Moose::Functor<Real> * const _material_density;
  const Moose::Functor<Real> * const _background_density;

  /// Current face volumetric-flux density.
  Real _face_volumetric_flux;

  /// Upwind interpolation coefficients for current face density.
  std::pair<Real, Real> _density_interp_coeffs;

  /// Optional material-compression mass-flux density on the current face.
  Real _face_explicit_correction;

  /// Frozen explicit correction used by the last assembled discrete mass equation.
  std::unordered_map<dof_id_type, Real> _assembled_explicit_correction;
};
