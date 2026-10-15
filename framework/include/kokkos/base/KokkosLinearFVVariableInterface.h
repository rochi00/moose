//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosSystem.h"

#include "MooseTypes.h"

#include <petscvec.h>

#include <set>
#include <vector>

class FEProblemBase;
class MooseObject;
class LinearFVGradientReader;

namespace Moose::Kokkos
{

/**
 * Device handle to the cell values of a linear finite volume variable in one solution state.
 * Values of the one-layer ghost cells are available.
 */
struct LinearFVCellValues
{
  /// Number of the variable's system
  unsigned int sys = libMesh::invalid_uint;
  /// Number of the variable in its system
  unsigned int var = libMesh::invalid_uint;
  /// Tag of the solution state
  TagID tag = Moose::INVALID_TAG_ID;

#ifdef MOOSE_KOKKOS_SCOPE
  /**
   * Get the value in a cell
   * @param systems The Kokkos systems, from SystemHolder::kokkosSystems()
   * @param elem The contiguous ID of an owned or ghost cell
   */
  KOKKOS_FUNCTION Real operator()(const Array<System> & systems,
                                  const ContiguousElementID elem) const
  {
    const auto & system = systems[sys];
    return system.getVectorDofValue(system.getElemLocalDofIndex(elem, 0, var), tag);
  }
#endif
};

/**
 * Device handle to the current cell gradients of a linear finite volume variable, as computed by
 * the variable's gradient method. The component arrays are set by
 * LinearFVVariableInterface::acquireGradients() before a dispatch. Gradients of the one-layer ghost
 * cells are available.
 */
struct LinearFVCellGradient
{
  /// Number of the variable's system
  unsigned int sys = libMesh::invalid_uint;
  /// Number of the variable in its system
  unsigned int var = libMesh::invalid_uint;
  /// Device arrays of the gradient components, indexed by local degree of freedom; nullptr past
  /// the mesh dimension
  const Real * components[3] = {nullptr, nullptr, nullptr};

#ifdef MOOSE_KOKKOS_SCOPE
  /**
   * Get the gradient in a cell
   * @param systems The Kokkos systems, from SystemHolder::kokkosSystems()
   * @param elem The contiguous ID of an owned or ghost cell
   */
  KOKKOS_FUNCTION Real3 operator()(const Array<System> & systems,
                                   const ContiguousElementID elem) const
  {
    const auto dof = systems[sys].getElemLocalDofIndex(elem, 0, var);
    Real3 gradient;
    for (unsigned int d = 0; d < 3; ++d)
      if (components[d])
        gradient(d) = components[d][dof];
    return gradient;
  }
#endif
};

/**
 * Interface for Kokkos linear finite volume objects that read the cell values of linear finite
 * volume variables, of any system and solution state. It records the systems an object reads, so
 * that the assembly synchronizes their solution vectors.
 */
class LinearFVVariableInterface
{
public:
  /**
   * Constructor
   * @param moose_object The object
   * @param problem The problem
   * @param own_system The number of the system of the object's own variable, which is always read
   */
  LinearFVVariableInterface(const MooseObject & moose_object,
                            FEProblemBase & problem,
                            const unsigned int own_system);
  /**
   * Copy constructor for parallel dispatch, which does not need the host bookkeeping
   */
  LinearFVVariableInterface(const LinearFVVariableInterface & object);

  /**
   * Get the systems whose solution vectors this object reads (host only)
   */
  const std::set<unsigned int> & readSystems() const { return _read_systems; }

  /**
   * Point the gradient handles of this object at the current gradients, before dispatching work
   * that reads them. releaseGradients() must follow once that work has finished.
   */
  void acquireGradients();
  /// Release the gradient arrays acquired by acquireGradients()
  void releaseGradients();

protected:
  /**
   * Get the cell values of a linear finite volume variable
   * @param name The variable name
   * @param state The solution state: 0 for the current solution, 1 for the old one, 2 for the
   * older one
   * @returns The device handle
   */
  LinearFVCellValues getLinearFVCellValues(const VariableName & name, const unsigned int state = 0);

  /**
   * Set up a handle to the current cell gradients of a linear finite volume variable. The handle
   * must be a member of this object, as acquireGradients() updates it.
   * @param name The variable name
   * @param gradient The handle
   * @param method_name The gradient method, or empty for the variable's default method
   */
  void requestLinearFVCellGradient(const VariableName & name,
                                   LinearFVCellGradient & gradient,
                                   const GradientMethodName & method_name = "");

private:
  /// The object
  const MooseObject & _lfv_moose_object;
  /// The problem
  FEProblemBase & _lfv_problem;
  /// Systems whose solution vectors this object reads (host only)
  std::set<unsigned int> _read_systems;

  /// A gradient handle of this object and the arrays it was pointed at
  struct GradientAccess
  {
    LinearFVCellGradient * handle;
    const LinearFVGradientReader * reader;
    std::vector<Vec> global;
    std::vector<Vec> local;
    std::vector<const PetscScalar *> arrays;
    /// Device copies of host arrays the device cannot address
    std::vector<Array<Real>> staged;
  };
  /// Gradient handles of this object (host only)
  std::vector<GradientAccess> _gradients;
};

} // namespace Moose::Kokkos
