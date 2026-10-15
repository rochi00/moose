//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

/**
 * Test-only postprocessor that checks the Kokkos finite volume mesh data. On host it compares every
 * face and cell value with the MOOSE FaceInfo and ElemInfo objects. On device it checks that the
 * faces of every owned cell close its surface, which exercises the cell-to-face table, the face
 * orientation and the device copies. Returns the number of failed checks over all processes.
 */
class KokkosFVMeshTest : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  KokkosFVMeshTest(const InputParameters & parameters);

  virtual void initialize() override { _failures = 0; }
  virtual void execute() override;
  virtual void finalize() override { gatherSum(_failures); }
  virtual Real getValue() const override { return _failures; }

private:
  /// Check the data against FaceInfo and ElemInfo on host
  void checkHost();
  /// Check the closure of the faces of every owned cell on device
  void checkDevice();

  /// Relative tolerance of the closure check, for the round-off of summing the area-weighted
  /// normals of a cell
  const Real _closure_tolerance;

  /// Number of failed checks
  Real _failures = 0;
};
