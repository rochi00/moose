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
  KOKKOS_FUNCTION Real faceConductance(const FaceDatum & datum, const Real diffusivity) const;

  /// Diffusion coefficient
  Moose::Kokkos::ReferenceWrapper<const KokkosParsedFunction> _diffusion_coeff;
  /// Whether the normal gradient on a face is corrected for the non-orthogonality of the mesh
  const bool _use_nonorthogonal_correction;
  /// The cell gradients of the variable, for the non-orthogonal correction
  Moose::Kokkos::LinearFVCellGradient _gradient;
};

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVDiffusion::InternalFaceContribution
KokkosLinearFVDiffusion::computeInternalFaceContribution(const FaceDatum & datum) const
{
  const auto diffusivity = _diffusion_coeff->value(_t, datum.faceCentroid());
  const auto conductance = faceConductance(datum, diffusivity);
  if (!_use_nonorthogonal_correction)
    return {conductance, -conductance, 0, 0};

  // Explicit part of the normal gradient: the face gradient, interpolated linearly from the cell
  // gradients, along the part of the unit normal that the implicit difference along the
  // cell-to-cell vector does not capture. The correction vector does not depend on the length of
  // the cell-to-cell vector, which therefore needs no normalization.
  const auto & systems = kokkosSystems();
  const auto elem = datum.elem();
  const auto neighbor = datum.neighbor();
  const auto normal = datum.faceNormal();
  const auto d_cn = datum.neighborCentroid() - datum.elemCentroid();
  const auto correction_vector = normal - 1 / (normal * d_cn) * d_cn;
  const auto gc = datum.faceGC();
  const auto face_gradient = gc * _gradient(systems, elem) + (1 - gc) * _gradient(systems, neighbor);
  const auto correction = diffusivity * (face_gradient * correction_vector) * datum.faceArea();
  return {conductance, -conductance, correction, -correction};
}

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVDiffusion::BoundaryFaceContribution
KokkosLinearFVDiffusion::computeBoundaryFaceContribution(
    const FaceDatum & datum, const int bc_index, const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const
{
  // A boundary condition prescribing a diffusive flux already includes the diffusion coefficient
  const auto diffusivity = _diffusion_coeff->value(_t, datum.faceCentroid());
  const auto face_flux_coefficient = boundaryIncludesMaterialPropertyMultiplier(bc_index)
                                         ? datum.faceArea()
                                         : diffusivity * datum.faceArea();
  auto rhs = face_flux_coefficient * boundaryNormalGradientSource(bc_index, bc_face_index);

  // The boundary conditions that compute the boundary value from the cell value, rather than
  // prescribing the flux, need the explicit part of the normal gradient as on internal faces, with
  // the cell gradient and the cell-to-face vector
  if (_use_nonorthogonal_correction && boundaryNeedsNonorthogonalCorrection(bc_index))
  {
    const auto normal = datum.faceNormal();
    const auto e_cf = datum.cellToFaceVector();
    const auto correction_vector = normal - 1 / (normal * e_cf) * e_cf;
    rhs += diffusivity * (_gradient(kokkosSystems(), datum.elem()) * correction_vector) *
           datum.faceArea();
  }

  return {-face_flux_coefficient * boundaryNormalGradientCoefficient(bc_index, bc_face_index), rhs};
}

KOKKOS_FUNCTION inline Real
KokkosLinearFVDiffusion::faceConductance(const FaceDatum & datum, const Real diffusivity) const
{
  // With the non-orthogonal correction, the implicit part of the normal gradient is the difference
  // of the cell values over the normal projection of the cell-to-cell distance
  if (!_use_nonorthogonal_correction)
    return diffusivity * datum.faceArea() / datum.faceDCNMag();

  const auto d_cn = datum.neighborCentroid() - datum.elemCentroid();
  const auto distance = d_cn * datum.faceNormal();
  return diffusivity / (distance < 0 ? -distance : distance) * datum.faceArea();
}
