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

class LinearSystem;

/**
 * An advection kernel that implements the advection term for the passive scalar transport equation.
 */
class LinearFVMultiPhaseFractionAdvection : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();
  LinearFVMultiPhaseFractionAdvection(const InputParameters & params);

  virtual Real computeElemMatrixContribution() override;

  virtual Real computeNeighborMatrixContribution() override;

  virtual Real computeElemRightHandSideContribution() override;

  virtual Real computeNeighborRightHandSideContribution() override;

  virtual Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;

  virtual Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;

  virtual void setupFaceData(const FaceInfo * face_info) override;

  /// Apply the semi-implicit MULES correction after the upwind solve.
  /// Called from the executioner after the linear system solve completes.
  /// @param system The linear system that was just solved (contains alpha_upwind)
  void applyMULESCorrection(LinearSystem & system);

  /// Whether semi-implicit MULES mode is active
  bool isSemiImplicitMULES() const { return _semi_implicit_mules; }

protected:
  /// Function to compute compression mass flux
  Real computeCompressionVelocityMassFlux();

  /// Function to compute the contribution from the compression velocity
  /// to the mass fluxes for the RHS term due to the non-orthogonal correction
  Real computeCompressionVelocityMassFluxNonOrthogonalRHSContribution();

  /// Function to get face value with low-order interpolation
  Real getLowOrderFaceValue(MooseLinearVariableFV<Real> & variable);

  /// Function to get face value with high-order interpolation
  Real getHighOrderFaceValue(MooseLinearVariableFV<Real> & variable);

  /// The Rhie-Chow user object that provides us with the face velocity
  RhieChowMassFluxMultiPhase & _mass_flux_provider;

  /// The dimension of the mesh
  const unsigned int _dim;

  /// The compression velocity control parameter
  const Real _c_alpha;

  /// The phase density used for compression velocity
  const Moose::Functor<Real> * _rho;

  /// Switch to enable/disable nonorthogonal correction in the stress term
  const bool _use_nonorthogonal_correction;

  /// Switch to activate MULES
  const bool _use_mules;

  /// Semi-implicit MULES mode: two-pass (implicit upwind + explicit limited correction)
  const bool _semi_implicit_mules;

  /// Use volumetric flux (phi) instead of mass flux (rho*phi) for advection
  const bool _use_volumetric_flux;

  /// Number of internal MULES iterations to perform
  const unsigned int _MULES_iterations;

  /// Enumerator with limiting methods
  enum LimiterMethod
  {
    MIN_MOD,
    VANLEER,
    VANALBADA,
    QUICK,
    VENKATAKRISHNAN,
    AVERAGE,
    UPWIND,
    UNKNOWN
  };

  // Function to map MooseEnum value to LimiterMethod
  LimiterMethod getLimiterMethod(const MooseEnum & limiter_enum)
  {
    std::string method = limiter_enum;
    if (method == "min_mod")
      return MIN_MOD;
    else if (method == "vanLeer")
      return VANLEER;
    else if (method == "vanAlbada")
      return VANALBADA;
    else if (method == "quick")
      return QUICK;
    else if (method == "venkatakrishnan")
      return VENKATAKRISHNAN;
    else if (method == "average")
      return AVERAGE;
    else if (method == "upwind")
      return UPWIND;
    else
      return UNKNOWN;
  }


private:
  /// Container for the current advected interpolation coefficients on the face to make sure
  /// we don't compute it multiple times for different terms.
  std::pair<Real, Real> _advected_interp_coeffs;

  /// Container for the velocity on the face which will be reused in the advection term's
  /// matrix and right hand side contribution
  Real _total_adv_mass_face_flux;

  /// The interpolation method to use for the advected quantity
  Moose::FV::InterpMethod _advected_interp_method;

  /// The limiter method
  MooseEnum _limiter_method;

  /// Face argument for higher order and low order face interpolation
  Moose::FaceArg _low_order_face;
  Moose::FaceArg _high_order_face;

  /// Container for the current advected interpolation coefficients for the compression velocity
  std::pair<Real, Real> _compression_interp_coeffs;

  /// Container for the compression velocity on the face
  Real _total_comp_mass_face_flux;

  /// CMULES face limiter
  Real _lambda_f;
};
