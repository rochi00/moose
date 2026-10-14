//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "WCNSLinearFVFlowPhysics.h"
#include "NS.h"

/**
 * Linear-FV segregated diffuse-interface flow physics for velocity-component momentum unknowns.
 *
 * This class reuses the common segregated Linear FV flow setup and overrides the diffuse-interface
 * pieces needed for large-density-ratio conservative coupling.
 */
class WCNSLinearFVConservativeDiffuseInterfaceFlowPhysics final : public WCNSLinearFVFlowPhysics
{
public:
  static InputParameters validParams();

  WCNSLinearFVConservativeDiffuseInterfaceFlowPhysics(const InputParameters & parameters);

private:
  /// The pressure equation of the reduced pressure: raw face rAUf and the explicit-force-free flux
  void addPressureCorrectionKernels() override;
  /// Pressure outlets with backflow handling, and the total-pressure condition on p_rgh
  void addOutletBC() override;
  /// No-slip walls carry a constrained zero-normal-flux pressure condition
  void addWallsBC() override;
  /// The diffuse-interface Rhie-Chow object, with gravity, reference point and face forces
  void addRhieChowUserObjects() override;
  void addFunctorMaterials() override;

  MooseFunctorName generatedGeometryFunctorName(const std::string & base_name) const;
  /// Creates the material publishing the interface curvature and capillary force density.
  void addSurfaceTensionFunctorMaterial();
  bool shouldCreateGeometryFunctorMaterial() const;

  unsigned short getNumberAlgebraicGhostingLayersNeeded() const override;

  const bool _create_geometry_functors;
};
