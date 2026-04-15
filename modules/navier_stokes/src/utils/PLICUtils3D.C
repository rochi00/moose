//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PLICUtils3D.h"

#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"

#include <cmath>
#include <algorithm>
#include <numeric>

namespace NS
{
namespace PLIC
{

namespace // anonymous helpers
{

/**
 * Compute area of a planar polygon using Newell's method.
 * Works for any 3D polygon; returns the signed area magnitude.
 */
Real
polygonArea(const std::vector<Point> & poly)
{
  if (poly.size() < 3)
    return 0.0;

  // Newell's method: area_vec = 0.5 * sum_i (v_i x v_{i+1})
  VectorValue<Real> cross(0, 0, 0);
  const auto n = poly.size();
  for (std::size_t i = 0; i < n; ++i)
  {
    const auto & vi = poly[i];
    const auto & vj = poly[(i + 1) % n];
    cross(0) += (vi(1) - vj(1)) * (vi(2) + vj(2));
    cross(1) += (vi(2) - vj(2)) * (vi(0) + vj(0));
    cross(2) += (vi(0) - vj(0)) * (vi(1) + vj(1));
  }
  return 0.5 * cross.norm();
}

/**
 * Clip a convex polygon by the half-space n . x <= d using
 * Sutherland-Hodgman. Returns the clipped polygon vertices.
 */
std::vector<Point>
clipPolygon(const std::vector<Point> & poly,
            const VectorValue<Real> & normal,
            Real offset)
{
  if (poly.empty())
    return {};

  std::vector<Point> output;
  output.reserve(poly.size() + 1);

  const auto n = poly.size();
  for (std::size_t i = 0; i < n; ++i)
  {
    const auto & curr = poly[i];
    const auto & next = poly[(i + 1) % n];

    const Real d_curr = normal * curr - offset; // positive = outside
    const Real d_next = normal * next - offset;

    if (d_curr <= 0.0) // curr is inside
    {
      output.push_back(curr);
      if (d_next > 0.0) // next is outside — add intersection
      {
        const Real t = d_curr / (d_curr - d_next);
        output.push_back(curr + t * (next - curr));
      }
    }
    else // curr is outside
    {
      if (d_next <= 0.0) // next is inside — add intersection
      {
        const Real t = d_curr / (d_curr - d_next);
        output.push_back(curr + t * (next - curr));
      }
    }
  }

  return output;
}

/**
 * Compute the volume of a convex polyhedron given as a set of outward-oriented
 * polygonal faces. Each face is fan-triangulated from its first vertex, and
 * each triangle contributes (1/6) * v0 . (v1 x v2) to the signed volume
 * (the standard divergence-theorem formula for triangulated surfaces).
 */
Real
polyhedronVolume(const std::vector<std::vector<Point>> & face_polys)
{
  Real vol = 0.0;
  for (const auto & poly : face_polys)
  {
    if (poly.size() < 3)
      continue;

    // Fan triangulation from vertex 0
    const auto & v0 = poly[0];
    for (std::size_t i = 1; i + 1 < poly.size(); ++i)
    {
      const auto & v1 = poly[i];
      const auto & v2 = poly[i + 1];
      // Signed volume of tetrahedron (origin, v0, v1, v2) = v0 . (v1 x v2) / 6
      vol += v0 * v1.cross(v2);
    }
  }
  return std::abs(vol) / 6.0;
}

/**
 * Brent's method to find root of f(x) = 0 in [a, b].
 * f(a) and f(b) must have opposite signs.
 */
Real
brentSolve(const std::function<Real(Real)> & f, Real a, Real b, Real tol, unsigned int max_iter)
{
  Real fa = f(a);
  Real fb = f(b);

  if (fa * fb > 0.0)
    return 0.5 * (a + b); // No bracket — return midpoint

  if (std::abs(fa) < tol)
    return a;
  if (std::abs(fb) < tol)
    return b;

  // Ensure |f(a)| >= |f(b)|
  if (std::abs(fa) < std::abs(fb))
  {
    std::swap(a, b);
    std::swap(fa, fb);
  }

  Real c = a, fc = fa;
  bool mflag = true;
  Real d = 0.0; // Previous step
  Real s;

  for (unsigned int iter = 0; iter < max_iter; ++iter)
  {
    if (std::abs(fb) < tol)
      return b;
    if (std::abs(b - a) < tol)
      return b;

    if (std::abs(fa - fc) > tol && std::abs(fb - fc) > tol)
    {
      // Inverse quadratic interpolation
      s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) +
          c * fa * fb / ((fc - fa) * (fc - fb));
    }
    else
    {
      // Secant method
      s = b - fb * (b - a) / (fb - fa);
    }

    // Conditions for bisection
    const Real s_min = (3.0 * a + b) / 4.0;
    const Real s_max = b;
    const bool cond1 =
        !((s > std::min(s_min, s_max)) && (s < std::max(s_min, s_max)));
    const bool cond2 = mflag && (std::abs(s - b) >= std::abs(b - c) / 2.0);
    const bool cond3 = !mflag && (std::abs(s - b) >= std::abs(c - d) / 2.0);
    const bool cond4 = mflag && (std::abs(b - c) < tol);
    const bool cond5 = !mflag && (std::abs(c - d) < tol);

    if (cond1 || cond2 || cond3 || cond4 || cond5)
    {
      s = (a + b) / 2.0;
      mflag = true;
    }
    else
      mflag = false;

    const Real fs = f(s);
    d = c;
    c = b;
    fc = fb;

    if (fa * fs < 0.0)
    {
      b = s;
      fb = fs;
    }
    else
    {
      a = s;
      fa = fs;
    }

    // Ensure |f(a)| >= |f(b)|
    if (std::abs(fa) < std::abs(fb))
    {
      std::swap(a, b);
      std::swap(fa, fb);
    }
  }

  return b;
}

} // anonymous namespace

// ============================================================================
// Public API
// ============================================================================

Real
clippedFaceArea(const std::vector<Point> & face_verts,
                const VectorValue<Real> & normal,
                Real offset)
{
  const auto clipped = clipPolygon(face_verts, normal, offset);
  return polygonArea(clipped);
}

Real
clippedVolume(const std::vector<Point> & vertices,
              const std::vector<std::vector<unsigned int>> & faces,
              const VectorValue<Real> & normal,
              Real offset)
{
  // Strategy: use the divergence theorem V = (1/6) * Σ v0·(v1×v2) over all
  // triangles of a closed surface. The clipped polyhedron is bounded by:
  //   (a) the clipped portions of the original faces, plus
  //   (b) the cap face on the clipping plane (n·x = offset).
  //
  // Rather than explicitly constructing the cap polygon (which requires
  // deduplication and angle-sorting), we use a coordinate trick:
  //
  // Choose a reference point P on the clipping plane: P = offset * n.
  // Then V = (1/6) |Σ (v0-P)·((v1-P)×(v2-P))| over ALL triangles of the
  // closed surface. For triangles on the cap face, every vertex satisfies
  // n·v = offset, so (v-P) is perpendicular to n for all cap vertices.
  // This means (v0-P)·((v1-P)×(v2-P)) = 0 for every cap triangle.
  //
  // Therefore: V = (1/6) |Σ_clipped_original_faces (v0-P)·((v1-P)×(v2-P))|
  // No cap construction needed!

  const Point ref_pt = offset * Point(normal(0), normal(1), normal(2));

  Real vol = 0.0;
  for (const auto & face : faces)
  {
    // Build face polygon from vertex indices
    std::vector<Point> face_poly;
    face_poly.reserve(face.size());
    for (auto idx : face)
      face_poly.push_back(vertices[idx]);

    // Clip face by the half-space n·x <= offset
    auto clipped = clipPolygon(face_poly, normal, offset);
    if (clipped.size() < 3)
      continue;

    // Fan-triangulate and accumulate signed volume relative to ref_pt
    const Point a0 = clipped[0] - ref_pt;
    for (std::size_t i = 1; i + 1 < clipped.size(); ++i)
    {
      const Point a1 = clipped[i] - ref_pt;
      const Point a2 = clipped[i + 1] - ref_pt;
      vol += a0 * a1.cross(a2);
    }
  }

  return std::abs(vol) / 6.0;
}

Real
computePlaneOffset3D(const std::vector<Point> & vertices,
                     const std::vector<std::vector<unsigned int>> & faces,
                     Real cell_volume,
                     const VectorValue<Real> & normal,
                     Real alpha,
                     Real tol)
{
  alpha = std::max(0.0, std::min(1.0, alpha));

  if (vertices.empty())
    return 0.0;

  // Compute signed projections h_i = n . x_i for all vertices
  std::vector<Real> projections(vertices.size());
  for (std::size_t i = 0; i < vertices.size(); ++i)
    projections[i] = normal * vertices[i];

  const Real h_min = *std::min_element(projections.begin(), projections.end());
  const Real h_max = *std::max_element(projections.begin(), projections.end());

  // Trivial cases
  if (alpha <= 0.0)
    return h_min;
  if (alpha >= 1.0)
    return h_max;

  const Real target_vol = alpha * cell_volume;

  // Sort unique projections to define bracketing intervals
  std::vector<Real> sorted_h = projections;
  std::sort(sorted_h.begin(), sorted_h.end());
  // Remove near-duplicates
  auto last = std::unique(sorted_h.begin(),
                          sorted_h.end(),
                          [](Real a, Real b) { return std::abs(a - b) < 1e-14; });
  sorted_h.erase(last, sorted_h.end());

  // Evaluate clipped volume at each breakpoint to find bracketing interval
  std::vector<Real> volumes(sorted_h.size());
  for (std::size_t i = 0; i < sorted_h.size(); ++i)
    volumes[i] = clippedVolume(vertices, faces, normal, sorted_h[i]);

  // Find the interval [h_i, h_{i+1}] containing target_vol
  std::size_t bracket_lo = 0;
  for (std::size_t i = 0; i + 1 < sorted_h.size(); ++i)
  {
    if (volumes[i + 1] >= target_vol - tol)
    {
      bracket_lo = i;
      break;
    }
  }

  const Real d_lo = sorted_h[bracket_lo];
  const Real d_hi = (bracket_lo + 1 < sorted_h.size()) ? sorted_h[bracket_lo + 1] : h_max;

  // Use Brent's method within the bracketing interval.
  // The analytical cubic inversion (Dai & Tong) is complex to implement robustly
  // for all cell topologies. Brent's method with the geometric clippedVolume evaluator
  // is reliable and converges in ~15-20 iterations to machine precision.
  // The clippedVolume calls are O(N_faces * N_verts_per_face) which is small (6-12).
  auto residual = [&](Real d)
  { return clippedVolume(vertices, faces, normal, d) - target_vol; };

  return brentSolve(residual, d_lo, d_hi, tol * cell_volume, 100);
}

Real
truncatedVolumeFraction3D(const std::vector<Point> & vertices,
                          const std::vector<std::vector<unsigned int>> & faces,
                          Real cell_volume,
                          const VectorValue<Real> & normal,
                          Real offset)
{
  if (cell_volume < 1e-30)
    return 0.0;

  const Real vol = clippedVolume(vertices, faces, normal, offset);
  return std::max(0.0, std::min(1.0, vol / cell_volume));
}

void
elemGeometry(const Elem & elem,
             std::vector<Point> & vertices,
             std::vector<std::vector<unsigned int>> & faces)
{
  vertices.clear();
  faces.clear();

  // Extract corner vertices (ignore mid-edge/mid-face nodes for higher-order)
  unsigned int n_corners = 0;

  switch (elem.type())
  {
    // 2D elements
    case TRI3:
    case TRI6:
    case TRI7:
      n_corners = 3;
      break;

    case QUAD4:
    case QUAD8:
    case QUAD9:
      n_corners = 4;
      break;

    // 3D elements
    case TET4:
    case TET10:
    case TET14:
      n_corners = 4;
      break;

    case HEX8:
    case HEX20:
    case HEX27:
      n_corners = 8;
      break;

    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      n_corners = 6;
      break;

    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      n_corners = 5;
      break;

    default:
      // Unsupported element type — extract all nodes as vertices
      n_corners = elem.n_vertices();
      break;
  }

  vertices.resize(n_corners);
  for (unsigned int i = 0; i < n_corners; ++i)
    vertices[i] = elem.point(i);

  // Define face connectivity based on element type.
  // libMesh node numbering conventions:
  //
  // TRI3:  0-1-2 (single face, 2D)
  // QUAD4: 0-1-2-3 (single face, 2D — treated as edges for 2D PLIC)
  //
  // TET4:  vertices 0,1,2,3
  //   face 0: 0-1-2 (outward)
  //   face 1: 0-3-1
  //   face 2: 1-3-2
  //   face 3: 0-2-3
  //
  // HEX8:  vertices 0-7
  //   face 0: 0-3-2-1 (bottom, -z)
  //   face 1: 0-1-5-4 (front, -y)
  //   face 2: 1-2-6-5 (right, +x)
  //   face 3: 2-3-7-6 (back, +y)
  //   face 4: 0-4-7-3 (left, -x)
  //   face 5: 4-5-6-7 (top, +z)
  //
  // PRISM6: vertices 0-5 (tri at 0-1-2, tri at 3-4-5)
  //   face 0: 0-2-1       (bottom tri)
  //   face 1: 3-4-5       (top tri)
  //   face 2: 0-1-4-3     (quad)
  //   face 3: 1-2-5-4     (quad)
  //   face 4: 0-3-5-2     (quad)
  //
  // PYRAMID5: vertices 0-4 (quad base 0-1-2-3, apex 4)
  //   face 0: 0-3-2-1     (base quad)
  //   face 1: 0-1-4       (tri)
  //   face 2: 1-2-4       (tri)
  //   face 3: 2-3-4       (tri)
  //   face 4: 0-4-3       (tri)

  if (elem.dim() == 2)
  {
    // For 2D elements, the "faces" are edges (needed for area computation
    // in the PLIC context). We also add top and bottom faces for the
    // extruded prism interpretation.
    if (n_corners == 3) // TRI3
    {
      faces = {{0, 1}, {1, 2}, {2, 0}};
    }
    else if (n_corners == 4) // QUAD4
    {
      faces = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    }
    return;
  }

  // 3D face connectivity
  switch (elem.type())
  {
    case TET4:
    case TET10:
    case TET14:
      faces = {{0, 1, 2}, {0, 3, 1}, {1, 3, 2}, {0, 2, 3}};
      break;

    case HEX8:
    case HEX20:
    case HEX27:
      faces = {{0, 3, 2, 1}, {0, 1, 5, 4}, {1, 2, 6, 5},
               {2, 3, 7, 6}, {0, 4, 7, 3}, {4, 5, 6, 7}};
      break;

    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      faces = {{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {0, 3, 5, 2}};
      break;

    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      faces = {{0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {0, 4, 3}};
      break;

    default:
      // For unknown 3D types, we cannot define faces.
      // The caller should handle this gracefully.
      break;
  }
}

Real
timeIntegratedSubmergedArea(const std::vector<Point> & face_verts,
                            const VectorValue<Real> & normal,
                            Real d0,
                            Real d_half,
                            Real d1,
                            Real dt)
{
  if (dt <= 0.0)
    return 0.0;

  // Simpson's 3-point quadrature: integral = (dt/6) * (f_0 + 4*f_half + f_1)
  const Real A0 = clippedFaceArea(face_verts, normal, d0);
  const Real A_half = clippedFaceArea(face_verts, normal, d_half);
  const Real A1 = clippedFaceArea(face_verts, normal, d1);

  return (dt / 6.0) * (A0 + 4.0 * A_half + A1);
}

} // namespace PLIC
} // namespace NS
