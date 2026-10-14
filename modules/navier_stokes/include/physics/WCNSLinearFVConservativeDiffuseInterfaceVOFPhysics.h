//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "WCNSLinearFVScalarTransportPhysics.h"

class WCNSLinearFVConservativeDiffuseInterfaceVOFPhysics final
  : public WCNSLinearFVScalarTransportPhysics
{
public:
  static InputParameters validParams();

  WCNSLinearFVConservativeDiffuseInterfaceVOFPhysics(const InputParameters & parameters);

protected:
  void addInitialConditions() override;
  void addMaterials() override;
  void addUserObjects() override;

  unsigned short getNumberAlgebraicGhostingLayersNeeded() const override { return 2; }

private:
  void addScalarAdvectionKernels() override;
  void addScalarDiffusionKernels() override {}
  void addScalarSourceKernels() override {}
  /// A volume fraction prescribed at an inlet is the inherited fixed-value condition; the
  /// diffuse-interface path needs nothing added to it. A wall carries no volume fraction
  /// condition, so that one stays empty.
  void addScalarInletBC() override { WCNSLinearFVScalarTransportPhysics::addScalarInletBC(); }
  void addScalarWallBC() override {}
  void addScalarOutletBC() override;

  void addMixtureFunctorMaterial(const std::string & object_suffix,
                                 const MooseFunctorName & property_name,
                                 const MooseFunctorName & phase_1_name,
                                 const MooseFunctorName & phase_2_name,
                                 bool limit_phase_fraction);
};
