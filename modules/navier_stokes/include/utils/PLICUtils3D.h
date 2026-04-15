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

namespace libMesh
{
class Elem;
}

namespace NS
{
namespace PLIC
{

/**
 * Clip a convex polygon by a half-plane (n_hat . x <= offset) and return
 * the area of the clipped sub-polygon. Uses Sutherland-Hodgman clipping.
 *
 * @param face_verts  Ordered vertices of the convex polygon
 * @param normal      Clipping plane normal
 * @param offset      Clipping plane offset (n . x = offset)
 * @return Area of the submerged (n . x <= offset) portion
 */
Real clippedFaceArea(const std::vector<Point> & face_verts,
                     const VectorValue<Real> & normal,
                     Real offset);

/**
 * Clip a convex polyhedron by a plane (n_hat . x <= offset) and return
 * the volume of the submerged region. The polyhedron is specified by its
 * vertices and face connectivity (each face is a list of vertex indices,
 * ordered outward).
 *
 * @param vertices  Polyhedron vertex coordinates
 * @param faces     Face connectivity (each face = ordered vertex indices)
 * @param normal    Clipping plane normal
 * @param offset    Clipping plane offset
 * @return Volume of the submerged region
 */
Real clippedVolume(const std::vector<Point> & vertices,
                   const std::vector<std::vector<unsigned int>> & faces,
                   const VectorValue<Real> & normal,
                   Real offset);

/**
 * Compute the PLIC plane offset for a general convex polyhedron using the
 * analytical algorithm of Dai & Tong (2019).
 *
 * Given interface normal n_hat and target volume fraction alpha, finds offset d
 * such that the volume of the polyhedron below the plane (n_hat . x <= d)
 * equals alpha * cell_volume.
 *
 * Falls back to Brent's method if the analytical solver encounters degenerate
 * geometry (coplanar vertices, zero-volume slices, etc.).
 *
 * @param vertices     Polyhedron vertex coordinates
 * @param faces        Face connectivity
 * @param cell_volume  Total cell volume
 * @param normal       Unit interface normal
 * @param alpha        Target volume fraction in [0, 1]
 * @param tol          Tolerance for Brent fallback (default 1e-12)
 * @return Global plane offset d
 */
Real computePlaneOffset3D(const std::vector<Point> & vertices,
                          const std::vector<std::vector<unsigned int>> & faces,
                          Real cell_volume,
                          const VectorValue<Real> & normal,
                          Real alpha,
                          Real tol = 1e-12);

/**
 * Forward problem: compute volume fraction for a given PLIC plane in a
 * general polyhedron (for verification/testing).
 *
 * @param vertices     Polyhedron vertex coordinates
 * @param faces        Face connectivity
 * @param cell_volume  Total cell volume
 * @param normal       Plane normal
 * @param offset       Plane offset
 * @return Volume fraction in [0, 1]
 */
Real truncatedVolumeFraction3D(const std::vector<Point> & vertices,
                               const std::vector<std::vector<unsigned int>> & faces,
                               Real cell_volume,
                               const VectorValue<Real> & normal,
                               Real offset);

/**
 * Extract vertex coordinates and face connectivity from a libMesh Elem.
 * Handles: QUAD4 (2D), TRI3 (2D), HEX8, TET4, PRISM6, PYRAMID5.
 * Higher-order elements (HEX20, HEX27, TET10) use only corner vertices.
 *
 * @param elem      The libMesh element
 * @param vertices  [out] Vertex coordinates
 * @param faces     [out] Face connectivity (each face = list of vertex indices,
 *                  ordered with outward-pointing normal via right-hand rule)
 */
void elemGeometry(const Elem & elem,
                  std::vector<Point> & vertices,
                  std::vector<std::vector<unsigned int>> & faces);

/**
 * Compute the time-integrated geometric volume flux of the liquid phase
 * through a polygonal face using Simpson's 3-point quadrature.
 *
 * The PLIC plane is advected linearly in time: d(t) = d_0 - u_n * t,
 * so the submerged area A_sub(t) is at most quadratic in t for convex faces.
 * Simpson's rule is therefore exact (zero truncation error).
 *
 * The returned flux is: (dt/6) * (A_0 + 4*A_half + A_1)
 * where A_i = clippedFaceArea(face_verts, normal, d_i).
 *
 * Multiply by the volumetric face flux Q_f to get the volume of liquid
 * transported through the face.
 *
 * @param face_verts  Ordered vertices of the face polygon
 * @param normal      PLIC plane normal
 * @param d0          Plane offset at t = 0
 * @param d_half      Plane offset at t = dt/2
 * @param d1          Plane offset at t = dt
 * @param dt          Timestep size
 * @return Time-integrated submerged area fraction (area * time)
 */
Real timeIntegratedSubmergedArea(const std::vector<Point> & face_verts,
                                 const VectorValue<Real> & normal,
                                 Real d0,
                                 Real d_half,
                                 Real d1,
                                 Real dt);

} // namespace PLIC
} // namespace NS
