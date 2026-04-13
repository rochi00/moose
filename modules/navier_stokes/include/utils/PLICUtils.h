//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"
#include "libmesh/vector_value.h"
#include "libmesh/point.h"

namespace NS
{
namespace PLIC
{

/**
 * Compute the PLIC plane offset for a 2D axis-aligned rectangle.
 *
 * Given a unit normal n_hat and volume fraction alpha in [0,1], finds the
 * signed offset d_local such that the volume below the line
 *   |n_x| * xi + |n_y| * eta = d_local
 * equals alpha * dx * dy, where (xi, eta) are local coordinates in [0,dx]x[0,dy].
 *
 * Uses the Scardovelli & Zaleski (2000) analytical inversion.
 *
 * @param n_hat  Unit interface normal
 * @param alpha  Volume fraction in [0, 1]
 * @param dx     Cell width in x
 * @param dy     Cell height in y
 * @return Local plane offset d_local
 */
Real computePlaneOffset2DRect(const VectorValue<Real> & n_hat, Real alpha, Real dx, Real dy);

/**
 * Forward problem: compute the volume fraction for a given PLIC plane in a 2D rectangle.
 * Used for verification.
 *
 * @param n_hat  Unit interface normal (absolute values used)
 * @param d_local Plane offset in local coordinates
 * @param dx     Cell width in x
 * @param dy     Cell height in y
 * @return Volume fraction alpha in [0, 1]
 */
Real truncatedVolumeFraction2DRect(const VectorValue<Real> & n_hat,
                                   Real d_local,
                                   Real dx,
                                   Real dy);

/**
 * Compute signed distance from a point to the PLIC plane n_hat . x = d.
 * Positive on the n_hat side (inside phase-1).
 */
Real signedDistanceToPlane(const Point & p, const VectorValue<Real> & n_hat, Real d);

/**
 * Convert a local PLIC offset (computed in the |n_x|, |n_y| frame on [0,dx]x[0,dy])
 * to a global offset for the plane n_hat . x = d_global, given the cell's lower-left corner.
 *
 * @param n_hat   Unit interface normal (original signs)
 * @param d_local Offset from computePlaneOffset2DRect
 * @param x0      Lower-left corner of the cell
 * @param dx      Cell width in x
 * @param dy      Cell height in y
 * @return Global offset d_global
 */
Real localToGlobalOffset(const VectorValue<Real> & n_hat,
                         Real d_local,
                         const Point & x0,
                         Real dx,
                         Real dy);

/**
 * Convert a global PLIC offset to local coordinates of a given rectangle.
 * Reverse of localToGlobalOffset.
 *
 * @param n_hat    Unit interface normal (original signs)
 * @param d_global Global offset for the plane n_hat . x = d_global
 * @param x0       Lower-left corner of the rectangle
 * @param dx       Rectangle width in x
 * @param dy       Rectangle height in y
 * @return Local offset d_local in the |n_x|, |n_y| frame on [0,dx]x[0,dy]
 */
Real globalToLocalOffset(const VectorValue<Real> & n_hat,
                         Real d_global,
                         const Point & x0,
                         Real dx,
                         Real dy);

/**
 * Compute the volume fraction of phase-1 in an arbitrary axis-aligned rectangle
 * using a PLIC plane defined in global coordinates.
 *
 * @param n_hat       Unit interface normal
 * @param d_global    Global plane offset (n_hat . x = d_global)
 * @param rect_origin Lower-left corner of the query rectangle
 * @param rect_dx     Width of the query rectangle
 * @param rect_dy     Height of the query rectangle
 * @return Volume fraction in [0, 1]
 */
Real volumeFractionInRect(const VectorValue<Real> & n_hat,
                          Real d_global,
                          const Point & rect_origin,
                          Real rect_dx,
                          Real rect_dy);

} // namespace PLIC
} // namespace NS
