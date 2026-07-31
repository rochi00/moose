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
 * Fixed-outer-iterate advection of specific enthalpy with the common face mass flux.
 *
 * The face enthalpy is reconstructed with the same bounded CUI procedure used for low-Mach mass
 * transport. Enthalpy is fixed at the current outer iterate, so this kernel contributes
 * -div(mass_flux * h_CUI) to the right hand side of the Newton temperature system.
 */
class LinearFVLowMachEnthalpyAdvection : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();

  LinearFVLowMachEnthalpyAdvection(const InputParameters & params);

  Real computeElemMatrixContribution() override;
  Real computeNeighborMatrixContribution() override;
  Real computeElemRightHandSideContribution() override;
  Real computeNeighborRightHandSideContribution() override;
  Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;
  Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;
  void setupFaceData(const FaceInfo * face_info) override;

  void beginFixedPointFluxCapture();
  void finishFixedPointFluxCapture();

  bool usesMassFlux(const std::string & mass_flux_name) const
  {
    return _mass_flux_functor && _mass_flux_functor_name == mass_flux_name &&
           !_mass_flux_is_integrated;
  }

  Real assembledFaceEnthalpyFlux(const FaceInfo & face_info) const;
  bool hasAssembledFaceEnthalpyFlux(const FaceInfo & face_info) const
  {
    return _assembled_face_contribution.count(face_info.id());
  }
  void clearAssembledFaceEnthalpyFluxes() { _assembled_face_contribution.clear(); }

protected:
  Real computeCUIFaceEnthalpy(const FaceInfo & face_info, bool upwind_is_elem) const;

  /// Specific enthalpy fixed at the current outer iterate.
  const Moose::Functor<Real> & _specific_enthalpy;

  /// Current dh/dT used to reconstruct the upwind enthalpy gradient.
  const Moose::Functor<Real> & _dh_dT;

  /// Prescribed specific enthalpy on boundary inflow.
  const Moose::Functor<Real> * const _inflow_specific_enthalpy;

  /// Optional Rhie-Chow provider of the common face mass flux.
  const RhieChowMassFlux * const _mass_flux_provider;

  /// Optional directly supplied common face mass-flux functor.
  const Moose::Functor<Real> * const _mass_flux_functor;

  /// Requested name of the directly supplied common face mass flux.
  const MooseFunctorName _mass_flux_functor_name;

  /// Whether the supplied mass flux already includes face area.
  const bool _mass_flux_is_integrated;

  /// Integrated mass-flux times specific-enthalpy contribution on the current face.
  Real _face_contribution;

  /// Integrated enthalpy flux fixed for the complete inner Newton solve.
  std::unordered_map<dof_id_type, Real> _fixed_point_face_contribution;

  /// Integrated enthalpy flux used by the last assembled Newton temperature equation.
  std::unordered_map<dof_id_type, Real> _assembled_face_contribution;

  bool _capture_fixed_point_flux = false;
  bool _use_fixed_point_flux = false;
};
