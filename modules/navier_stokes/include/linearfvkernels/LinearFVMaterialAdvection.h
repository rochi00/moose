//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "LinearFVFluxKernel.h"

class RhieChowMassFlux;

/**
 * Implicit finite-volume discretization of material advection by a generally non-divergence-free
 * velocity:
 *
 *   div(u q) - q div(u).
 *
 * Both terms are assembled face by face from the same volumetric flux. This makes the discrete
 * operator annihilate a constant q exactly, including when div(u) is nonzero.
 */
class LinearFVMaterialAdvection : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();

  LinearFVMaterialAdvection(const InputParameters & params);

  bool usesDeferredCorrection() const;
  const MooseFunctorName & volumetricFaceFluxName() const { return _volumetric_face_flux_name; }

  void addMatrixContribution() override;
  void addRightHandSideContribution() override;

  Real computeElemMatrixContribution() override;
  Real computeNeighborMatrixContribution() override;
  Real computeElemRightHandSideContribution() override;
  Real computeNeighborRightHandSideContribution() override;
  Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;
  Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;
  void setupFaceData(const FaceInfo * face_info) override;

protected:
  /// Volumetric face flux oriented with the FaceInfo normal.
  Real volumetricFaceFlux(const FaceInfo & face_info) const;

  /// Legacy bounded cubic-upwind face value for a uniform-grid internal face.
  Real computeCUIFaceValue(const FaceInfo & face_info, bool upwind_is_elem) const;

  /// Geometry-aware MUSCL value reconstructed at the actual subface centroid.
  Real computeMUSCLFaceValue(const FaceInfo & face_info, bool upwind_is_elem) const;

  /// Name of the directly supplied volumetric face-flux functor.
  const MooseFunctorName _volumetric_face_flux_name;

  /// Optional Rhie-Chow provider of the volumetric face flux.
  const RhieChowMassFlux * const _volumetric_flux_provider;

  /// Optional directly supplied volumetric face-flux functor.
  const Moose::Functor<Real> * const _volumetric_face_flux;

  /// Optional velocity components used to construct the oriented face flux.
  const Moose::Functor<Real> * const _u;
  const Moose::Functor<Real> * const _v;
  const Moose::Functor<Real> * const _w;

  /// Interpolation method for the advected field.
  Moose::FV::InterpMethod _advected_interp_method;

  /// Face interpolation coefficients multiplying the element and neighbor values.
  std::pair<Real, Real> _advected_interp_coeffs;

  /// Volumetric face-flux density oriented with the FaceInfo normal.
  Real _face_volumetric_flux;

  /// High-order minus donor-cell material-advection flux density.
  Real _face_deferred_correction;
};
