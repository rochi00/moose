//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosDispatcher.h"
#include "KokkosFVDatum.h"
#include "KokkosMesh.h"
#include "KokkosSystem.h"
#include "KokkosVariable.h"

#include "BoundaryRestrictableRequired.h"
#include "FunctionInterface.h"
#include "MeshChangedInterface.h"
#include "MooseObject.h"
#include "PostprocessorInterface.h"
#include "Restartable.h"
#include "SetupInterface.h"
#include "TransientInterface.h"
#include "UserObjectInterface.h"
#include "VectorPostprocessorInterface.h"

#include <map>
#include <optional>

class FEProblemBase;
class InputParameters;

namespace Moose::Kokkos
{

/**
 * Base class for Kokkos linear finite volume boundary conditions. Boundary conditions provide
 * boundary data relations for flux kernels; they do not assemble directly into tagged vectors or
 * matrices.
 */
class LinearFVBoundaryCondition : public MooseObject,
                                  public SetupInterface,
                                  public FunctionInterface,
                                  public UserObjectInterface,
                                  public TransientInterface,
                                  public PostprocessorInterface,
                                  public VectorPostprocessorInterface,
                                  public Restartable,
                                  public MeshChangedInterface,
                                  public BoundaryRestrictableRequired,
                                  public MeshHolder,
                                  public FVMeshHolder,
                                  public SystemHolder
{
public:
  static InputParameters validParams();

  LinearFVBoundaryCondition(const InputParameters & parameters);
  LinearFVBoundaryCondition(const LinearFVBoundaryCondition & object);

  /// Tag dispatch type for boundary value relation computation
  struct BoundaryValueLoop
  {
  };
  /// Tag dispatch type for boundary normal gradient relation computation
  struct BoundaryNormalGradientLoop
  {
  };

  /**
   * Affine boundary relation used by Kokkos linear FV kernels:
   *   boundary_quantity = coefficient * cell_unknown + source
   */
  struct BoundaryRelation
  {
    Real coefficient = 0;
    Real source = 0;
  };

  /// Boundary relation coefficient/source arrays indexed by BC-local boundary face
  struct BoundaryRelationData
  {
    Array<Real> coefficient;
    Array<Real> source;
  };

  /// Boundary relation data indexed by relation type and BC-local boundary face
  struct BoundaryData
  {
    BoundaryRelationData value;
    BoundaryRelationData normal_gradient;
  };

  virtual void initialSetup() override;

  /**
   * Get the Kokkos variable this boundary condition supplies data for
   * @returns The Kokkos variable
   */
  Variable variable() const { return _var; }

  /**
   * Build the faces of this boundary condition and its relation arrays, if not built yet. The
   * faces are those of the face ranges of the Kokkos finite volume mesh where the variable lives on
   * one side only, the cell on that side is owned, and the single boundary ID is one of this
   * boundary condition's. Flux kernels call it, as their setup can run before this object's.
   */
  void setupFaces();

  /**
   * Get the boundary data relation arrays owned by this boundary condition, indexed by face of
   * this boundary condition
   */
  const BoundaryData & boundaryData() const { return _boundary_data; }

  /// Number of faces of this boundary condition
  MOOSE_KOKKOS_INDEX_TYPE numBoundaryFaces() const { return _faces.size(); }

  /**
   * Get the index in the relation arrays of the first face of a face range of the Kokkos finite
   * volume mesh
   * @param range The index of the range in FVMesh::getFaceRanges()
   * @returns The index, or nothing if the range is not one of this boundary condition's
   */
  std::optional<MOOSE_KOKKOS_INDEX_TYPE> faceRangeOffset(const unsigned int range) const;

  /// Whether this boundary condition overrides the boundary value relation hook
  bool hasBoundaryValue() const;
  /// Whether this boundary condition overrides the boundary normal gradient relation hook
  bool hasBoundaryNormalGradient() const;

  /// Dispatch boundary value relation computation, without waiting for it to finish
  void computeBoundaryValueData();
  /// Dispatch boundary normal gradient relation computation, without waiting for it to finish
  void computeBoundaryNormalGradientData();

  /**
   * Default boundary value relation hook. Derived boundary conditions should override this when
   * they can supply boundary values.
   */
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryValue(const FaceDatum &) const
  {
    ::Kokkos::abort("Default computeBoundaryValue() should never be called. Make sure you "
                    "properly redefined this method in your class without typos.");

    return {};
  }

  /**
   * Default boundary normal gradient relation hook. Derived boundary conditions should override
   * this when they can supply boundary normal gradients.
   */
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryRelation computeBoundaryNormalGradient(const FaceDatum &) const
  {
    ::Kokkos::abort("Default computeBoundaryNormalGradient() should never be called. Make sure "
                    "you properly redefined this method in your class without typos.");

    return {};
  }

  /**
   * Functions used to check whether derived boundary conditions override the boundary data hooks.
   */
  ///@{
  template <typename Derived>
  static auto defaultBoundaryValue()
  {
    return &LinearFVBoundaryCondition::computeBoundaryValue<Derived>;
  }

  template <typename Derived>
  static auto defaultBoundaryNormalGradient()
  {
    return &LinearFVBoundaryCondition::computeBoundaryNormalGradient<Derived>;
  }
  ///@}

  /**
   * Boundary value dispatch loop body; writes the per-face boundary value relation
   * @param tid The thread ID of the current boundary side
   * @param bc The concrete boundary condition object
   */
  template <typename Derived>
  KOKKOS_FUNCTION void operator()(BoundaryValueLoop, const ThreadID tid, const Derived & bc) const;

  /**
   * Boundary normal gradient dispatch loop body; writes the per-face boundary normal gradient
   * relation
   * @param tid The thread ID of the current boundary side
   * @param bc The concrete boundary condition object
   */
  template <typename Derived>
  KOKKOS_FUNCTION void
  operator()(BoundaryNormalGradientLoop, const ThreadID tid, const Derived & bc) const;

protected:
  /// The face of a thread, seen from the cell where the variable lives
  KOKKOS_FUNCTION FaceDatum boundaryFaceDatum(const ThreadID tid) const
  {
    const auto range = _faces.range(tid);
    return FaceDatum(_faces.face(tid, range), _face_range_flips[range], kokkosFVMesh());
  }

  /// Reference to the finite element problem
  FEProblemBase & _fe_problem;

  /// Kokkos variable
  Variable _var;

  /// Boundary relation arrays owned by this boundary condition
  BoundaryData _boundary_data;

  /// Faces of this boundary condition
  FaceRangeList _faces;
  /// Whether the variable lives on the neighbor side of the faces of each range of _faces
  Array<bool> _face_range_flips;
  /// Index in _faces of each face range of the Kokkos finite volume mesh this boundary condition
  /// covers (host only)
  std::map<unsigned int, unsigned int> _mesh_range_to_range;
  /// Whether setupFaces() ran
  bool _faces_built = false;

  /// Dispatcher for boundary value computation
  std::unique_ptr<DispatcherBase> _boundary_value_dispatcher;
  /// Dispatcher for boundary normal gradient computation
  std::unique_ptr<DispatcherBase> _boundary_normal_gradient_dispatcher;

  /**
   * TODO: Move to TransientInterface
   */
  ///@{
  /// Current time
  Scalar<Real> _t;
  /// Old (previous time step) time
  Scalar<const Real> _t_old;
  /// Current time step number
  Scalar<int> _t_step;
  /// Current time step size
  Scalar<Real> _dt;
  /// Previous time step size
  Scalar<Real> _dt_old;
  ///@}
};

template <typename Derived>
KOKKOS_FUNCTION void
LinearFVBoundaryCondition::operator()(BoundaryValueLoop,
                                      const ThreadID tid,
                                      const Derived & bc) const
{
  const auto datum = boundaryFaceDatum(tid);
  const auto relation = bc.template computeBoundaryValue<Derived>(datum);
  _boundary_data.value.coefficient[tid] = relation.coefficient;
  _boundary_data.value.source[tid] = relation.source;
}

template <typename Derived>
KOKKOS_FUNCTION void
LinearFVBoundaryCondition::operator()(BoundaryNormalGradientLoop,
                                      const ThreadID tid,
                                      const Derived & bc) const
{
  const auto datum = boundaryFaceDatum(tid);
  const auto relation = bc.template computeBoundaryNormalGradient<Derived>(datum);
  _boundary_data.normal_gradient.coefficient[tid] = relation.coefficient;
  _boundary_data.normal_gradient.source[tid] = relation.source;
}

} // namespace Moose::Kokkos
