//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"

/**
 * Writes the finite volume cell and face data (ElemInfo/FaceInfo geometry and the dof index of
 * one linear FV variable) to a binary file. Used by the standalone Kokkos linear FV layout
 * benchmark so that it operates on exactly the mesh data MOOSE builds. Serial only.
 */
class FVMeshDumper : public GeneralUserObject
{
public:
  static InputParameters validParams();

  FVMeshDumper(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}
};
