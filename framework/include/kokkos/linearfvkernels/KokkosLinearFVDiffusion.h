//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosLinearFVFluxKernel.h"
#include "KokkosParsedFunction.h"

/**
 * Kokkos linear finite volume flux kernel implementing the diffusion term. It is the Kokkos analog
 * of LinearFVDiffusion.
 */
class KokkosLinearFVDiffusion : public Moose::Kokkos::LinearFVFluxKernel
{
public:
  static InputParameters validParams();

  KokkosLinearFVDiffusion(const InputParameters & parameters);

  virtual bool needsBoundaryNormalGradientData() const override { return true; }

  template <typename Derived>
  KOKKOS_FUNCTION InternalFaceContribution
  computeInternalFaceContribution(const FaceDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryFaceContribution computeBoundaryFaceContribution(
      const FaceDatum & datum, const int bc_index, const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const;

private:
  /**
   * Method for evaluating the face diffusion coefficient times the face area divided by the
   * distance between cell centers. The name "conductance" is used here for this finite volume face
   * coefficient because heat conduction problems provide a thermal conductivity as the diffusion
   * coefficient, which has conductance units after multiplying by area over length. Similar
   * "diffusive conductance" terminology is also used in mass transport contexts; see
   * https://www.goldsim.com/Courses/ContaminantTransport/Unit8/Lesson3/
   */
  KOKKOS_FUNCTION Real faceConductance(const FaceDatum & datum) const;

  /// Diffusion coefficient
  Moose::Kokkos::ReferenceWrapper<const KokkosParsedFunction> _diffusion_coeff;
};

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVDiffusion::InternalFaceContribution
KokkosLinearFVDiffusion::computeInternalFaceContribution(const FaceDatum & datum) const
{
  const auto conductance = faceConductance(datum);
  return {conductance, -conductance, 0, 0};
}

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVDiffusion::BoundaryFaceContribution
KokkosLinearFVDiffusion::computeBoundaryFaceContribution(
    const FaceDatum & datum, const int bc_index, const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const
{
  const auto face_flux_coefficient =
      _diffusion_coeff->value(_t, datum.faceCentroid()) * datum.faceArea();
  return {-face_flux_coefficient * boundaryNormalGradientCoefficient(bc_index, bc_face_index),
          face_flux_coefficient * boundaryNormalGradientSource(bc_index, bc_face_index)};
}

KOKKOS_FUNCTION inline Real
KokkosLinearFVDiffusion::faceConductance(const FaceDatum & datum) const
{
  return _diffusion_coeff->value(_t, datum.faceCentroid()) * datum.faceArea() / datum.faceDCNMag();
}
