//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosArray.h"
#include "KokkosFVMesh.h"

#include "MooseTypes.h"

#include <memory>
#include <string>
#include <vector>

class FEProblemBase;
class FVFaceInterpolationMethod;
class KokkosParsedFunction;
class MooseObject;

namespace Moose::Kokkos
{

class LinearFVVariableInterface;

/**
 * Device handle to a property of a Kokkos linear finite volume object: a constant, or its values on
 * the cells and faces the object needs. The arrays are set and filled by
 * LinearFVPropertyInterface::updateLinearFVProperties().
 */
struct LinearFVProperty
{
  /// The value of a constant property
  Real constant = 0;
  /// Values on the owned and ghost cells by contiguous element ID, or nullptr if constant
  const Real * cell_values = nullptr;
  /// Values on the faces of the object by contiguous face ID, or nullptr if constant
  const Real * face_values = nullptr;

#ifdef MOOSE_KOKKOS_SCOPE
  /// Get the value on a cell
  KOKKOS_FUNCTION Real cell(const ContiguousElementID elem) const
  {
    return cell_values ? cell_values[elem] : constant;
  }
  /// Get the value on a face
  KOKKOS_FUNCTION Real face(const ContiguousFaceID face) const
  {
    return face_values ? face_values[face] : constant;
  }
#endif
};

/**
 * Interface for Kokkos linear finite volume objects with properties given as functors: a number,
 * the name of a KokkosConstantFunction or KokkosParsedFunction, the name of a linear finite volume
 * variable, or the property_name of a KokkosParsedFunctorMaterial. As for the host functors,
 *
 * - a cell value is evaluated at the cell centroid, which is the vertex average of the element as
 *   ElemInfo::centroid() and Moose::ElemArg::getPoint(); a variable gives its cell value;
 * - a face value is evaluated at the face centroid, or interpolated from the cell values with an
 *   FVFaceInterpolationMethod (FVGeometricAverage or FVHarmonicAverage) on the faces with the
 *   object's variable on both sides. A function evaluated from the object's variable's side, as
 *   with the face argument of singleSidedFaceArg(), moves 1e-8 along the normal into the
 *   variable's cell on faces between blocks; with an unsided face argument, as functorFaceArg()
 *   gives for a function, it does not. A variable gives its face value: the linear interpolation
 *   of its cell values where it lives on both sides, else the boundary value of its boundary
 *   condition, else its cell value (MooseLinearVariableFV::evaluate()); an expression combines the
 *   face values of its inputs.
 *
 * Values of functions are computed on the device when the time changes, values of variables and
 * expressions of variables at every update, rather than in every face and cell pass.
 */
class LinearFVPropertyInterface
{
public:
  /**
   * Constructor
   * @param moose_object The object
   * @param problem The problem
   * @param variables The variable interface of the object, which reads the variables of the
   * properties
   */
  LinearFVPropertyInterface(const MooseObject & moose_object,
                            FEProblemBase & problem,
                            LinearFVVariableInterface & variables);
  /**
   * Copy constructor for parallel dispatch, which does not need the host bookkeeping
   */
  LinearFVPropertyInterface(const LinearFVPropertyInterface & object);
  ~LinearFVPropertyInterface();

protected:
  /// Where an object needs the values of a property
  struct PropertyUse
  {
    /// On the owned cells
    bool cells = false;
    /// On the faces with the variable on both sides
    bool internal_faces = false;
    /// On the faces with the variable on one side
    bool boundary_faces = false;
    /// Whether the functor is evaluated on the faces with the variable on one side without a side,
    /// rather than from the variable's side
    bool unsided_boundary_faces = false;
    /// Interpolation of the cell values to the faces with the variable on both sides, or nullptr
    /// to evaluate the functor on them
    const FVFaceInterpolationMethod * interpolation = nullptr;
  };

  /**
   * Set up a handle to a property. The handle must be a member of this object, as
   * updateLinearFVProperties() updates it. Numbers and functions are resolved here, variables
   * and functor material properties by setupLinearFVProperties().
   * @param param The parameter holding the functor name
   * @param property The handle
   * @param use Where the object needs the values
   */
  void requestLinearFVProperty(const std::string & param,
                               LinearFVProperty & property,
                               const PropertyUse & use);

  /**
   * Resolve the properties given by variables and functor material properties, once all objects
   * exist; called from initialSetup()
   */
  void setupLinearFVProperties();

  /**
   * Fill the property values, unless they are up to date. The work is dispatched without waiting
   * for it to finish, unless a property reads boundary values of a variable, whose boundary
   * conditions are updated first.
   * @param t The time
   * @param internal_faces The faces of this object with the variable on both sides, or nullptr
   * @param boundary_faces The faces of this object with the variable on one side, or nullptr
   * @param boundary_flips Whether the variable lives on the neighbor side of the faces of each
   * range of boundary_faces
   */
  void updateLinearFVProperties(Real t,
                                const FaceRangeList * internal_faces,
                                const FaceRangeList * boundary_faces,
                                const std::vector<bool> & boundary_flips);

private:
  /// The inputs, boundary conditions and face value sources of a property given by variables
  struct ExpressionData;

  /// A property of this object and its values
  struct PropertyData
  {
    /// The parameter holding the functor name
    std::string param;
    /// The handle
    LinearFVProperty * handle;
    /// Where the object needs the values
    PropertyUse use;
    /// The function, or nullptr
    const KokkosParsedFunction * function = nullptr;
    /// The variables and expression, or nullptr
    std::shared_ptr<ExpressionData> expression;
    /// Whether the property waits for setupLinearFVProperties()
    bool deferred = false;
    /// Values on the owned and ghost cells
    Array<Real> cell_values;
    /// Values on the faces
    Array<Real> face_values;
    /// Direction of the offset of the evaluation point on each range of boundary faces, kept with
    /// the values so that it outlives the fill
    Array<int> offset_sign;
    /// The time of the values
    Real time = 0;
    /// Whether the values were filled
    bool filled = false;
  };

  /**
   * Set up a property given by a linear finite volume variable or a KokkosParsedFunctorMaterial
   * @param data The property
   * @param name The functor name
   */
  void setupExpression(PropertyData & data, const std::string & name);

  /// The object
  const MooseObject & _lfvp_moose_object;
  /// The problem
  FEProblemBase & _lfvp_problem;
  /// The variable interface of the object (host only)
  LinearFVVariableInterface * _lfvp_variables;
  /// Properties of this object (host only)
  std::vector<PropertyData> _properties;
};

} // namespace Moose::Kokkos
