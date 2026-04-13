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
#include "MooseLinearVariableFV.h"

class LinearSystem;

namespace libMesh
{
class LinearImplicitSystem;
}

/**
 * Computes a PLIC (Piecewise Linear Interface Calculation) reconstruction
 * for each interfacial cell. The interface is represented as a plane
 *   n_hat . x = d
 * where n_hat comes from grad(phi) or grad(alpha) and d is determined
 * by matching the cell's volume fraction alpha.
 *
 * If level_set_variable is provided, normals come from grad(phi).
 * Otherwise, normals come from grad(alpha) (Youngs method).
 *
 * Currently supports 2D axis-aligned QUAD4 elements only.
 */
class PLICReconstruction : public GeneralUserObject
{
public:
  static InputParameters validParams();
  PLICReconstruction(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}

  struct PLICPlane
  {
    VectorValue<Real> n_hat;
    Real d;
  };

  /// Whether the element has a PLIC plane (is interfacial)
  bool hasPlane(dof_id_type elem_id) const;

  /// Get the PLIC plane for an element (must call hasPlane first)
  const PLICPlane & getPlane(dof_id_type elem_id) const;

  /// Get all reconstructed planes (for iteration)
  const std::unordered_map<dof_id_type, PLICPlane> & planes() const { return _planes; }

protected:
  /// The VOF variable (alpha)
  MooseLinearVariableFV<Real> & _alpha_var;

  /// Linear system for alpha (to compute gradients when no level-set)
  LinearSystem & _alpha_system;

  /// Whether a level-set variable is used for normals
  const bool _use_level_set;

  /// The level-set variable (phi) — gradient provides interface normal (optional)
  MooseLinearVariableFV<Real> * _phi_var;

  /// Linear system for phi (to compute gradients)
  LinearSystem * _phi_system;

  /// Tolerance for identifying interfacial cells
  const Real _alpha_tol;

  /// Per-cell PLIC plane data
  std::unordered_map<dof_id_type, PLICPlane> _planes;
};
