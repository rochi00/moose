//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MathFVUtils.h"
#include "NSFVUtils.h"
#include "MooseMesh.h"
#include "MooseFunctor.h"
#include "FaceInfo.h"
#include "ElemInfo.h"

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

namespace NS
{
/**
 * The FaceInfo covering a side of an element, or null when that side has none, which happens on a
 * remote side and on a side with neither a FaceInfo nor a neighbor.
 */
const FaceInfo *
faceInfoForElementSide(const MooseMesh & mesh, const Elem * elem, const unsigned int side);

/**
 * A functor's value on a face: a geometrically weighted average of the two cell values across an
 * internal face, and a one-sided evaluation on a boundary face.
 * @param side_elem the element whose side the face is, used for the one-sided evaluation
 * @param face_is_internal whether the face is internal to the domain the functor is defined on
 */
Real faceScalarValue(const Moose::Functor<Real> & functor,
                     const FaceInfo & fi,
                     const Elem * const side_elem,
                     const bool face_is_internal,
                     const Moose::StateArg & time_arg);

/**
 * A functor's face normal gradient, carrying the explicit non-orthogonal correction across an
 * internal face. On a boundary face the correction is zero and the boundary patch gradient is
 * what the face sees.
 * @param boundary_elem the element to evaluate the boundary patch gradient on
 * @param face_is_internal whether the face is internal to the domain the functor is defined on
 */
Real faceNormalScalarGradient(const Moose::Functor<Real> & functor,
                              const FaceInfo & fi,
                              const Elem * const boundary_elem,
                              const bool face_is_internal,
                              const Moose::StateArg & time_arg);

/**
 * The explicit non-orthogonal correction vector of an internal face: zero on an orthogonal face,
 * where the cell gradients it multiplies need not be evaluated.
 */
Point faceNonOrthogonalCorrection(const FaceInfo & fi);

/**
 * Face-normal gradient on an internal face from the two cell values and the already evaluated
 * cell gradients, so that a caller sweeping all faces evaluates each cell gradient once.
 */
Real faceNormalScalarGradient(const Moose::Functor<Real> & functor,
                              const FaceInfo & fi,
                              const Moose::StateArg & time_arg,
                              const RealVectorValue & elem_grad,
                              const RealVectorValue & neighbor_grad);

/**
 * A cell gradient fitted by least squares to the jumps between the cell value and its face values,
 * over the offsets from the cell centroid to the face centroids. Unlike Green-Gauss this does not
 * need a closed cell, so it degrades gracefully on a partially coupled patch: it returns a zero
 * gradient when fewer than \p dim usable sides remain.
 *
 * The predicate is a template parameter rather than a std::function because this runs per cell per
 * side.
 *
 * @param is_internal_face called with a FaceInfo, returns whether that face is internal to the
 *        domain the functor is defined on
 */
template <typename InternalFacePredicate>
RealVectorValue
leastSquaresCellGradient(const Moose::Functor<Real> & functor,
                         const ElemInfo & elem_info,
                         const MooseMesh & mesh,
                         const unsigned int dim,
                         const InternalFacePredicate & is_internal_face,
                         const Moose::StateArg & time_arg)
{
  const Elem * const elem = elem_info.elem();
  if (!elem)
    return RealVectorValue();

  const Real cell_value = functor(Moose::ElemArg{elem, false}, time_arg);
  std::vector<Point> offsets;
  std::vector<Real> increments;
  offsets.reserve(elem->n_sides());
  increments.reserve(elem->n_sides());

  for (const auto side : make_range(elem->n_sides()))
  {
    const FaceInfo * const fi_loc = faceInfoForElementSide(mesh, elem, side);
    if (!fi_loc)
      continue;

    const Point offset = fi_loc->faceCentroid() - elem_info.centroid();
    if (offset.norm() <= libMesh::TOLERANCE)
      continue;

    offsets.push_back(offset);
    increments.push_back(
        faceScalarValue(functor, *fi_loc, elem, is_internal_face(*fi_loc), time_arg) - cell_value);
  }

  if (offsets.size() < dim)
    return RealVectorValue();

  DenseMatrix<Real> matrix(offsets.size(), dim);
  DenseVector<Real> rhs(offsets.size());
  DenseVector<Real> solution(dim);

  for (const auto row : index_range(offsets))
  {
    rhs(row) = increments[row];
    for (const auto dim_i : make_range(dim))
      matrix(row, dim_i) = offsets[row](dim_i);
  }

  matrix.svd_solve(rhs, solution);

  RealVectorValue gradient;
  for (const auto dim_i : make_range(dim))
    gradient(dim_i) = solution(dim_i);

  return gradient;
}
}
