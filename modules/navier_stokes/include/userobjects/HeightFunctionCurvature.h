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
#include <unordered_map>

class PLICReconstruction;
class LinearSystem;

namespace libMesh
{
class LinearImplicitSystem;
}

/**
 * 2D weighted least-squares parabola fit in a local coordinate frame.
 * Fits eta = a[0]*xi^2 + a[1]*xi + a[2] where (xi, eta) are local coords
 * aligned with the interface (xi tangent, eta normal).
 * Based on Basilisk's parabola.h (Popinet 2009).
 */
struct ParabolaFit2D
{
  /// Origin of local frame (interface centroid)
  Point origin;
  /// Local frame: tangent (t) and normal (n) directions
  Real tx, ty, nx, ny;
  /// Normal equations matrix (symmetric 3x3)
  Real M[3][3];
  /// Right-hand side
  Real rhs[3];
  /// Fitted coefficients
  Real a[3];

  /// Initialize with origin and normal direction
  void init(const Point & o, const VectorValue<Real> & normal);
  /// Add a weighted point
  void addPoint(const Point & p, Real weight);
  /// Solve the normal equations. Returns false if singular.
  bool solve();
  /// Compute curvature from the fit. Clamps to |kappa| <= kappa_max.
  Real curvature(Real kappa_max) const;
};

/**
 * Computes interface curvature using height functions (Popinet 2009).
 *
 * Height values match Basilisk's convention: signed distance from cell
 * center to interface, normalized by Delta (dimensionless). Orientation
 * is encoded via HSHIFT = 20.
 *
 * Algorithm:
 *   1. Compute heights for all cells (two-pass half-column, both directions)
 *   2. Column propagation (extend heights to neighbors)
 *   3. Curvature from finite differences: kappa = h_xx / (1 + h_x^2)^{3/2}
 *      with sign flip based on c[+1] - c[-1]
 *
 * Fallback hierarchy (matching Basilisk):
 *   1. Standard height-function curvature (try both directions)
 *   2. Neighbor-averaged curvature from valid HF neighbors
 *   3. Parabolic fit of PLIC interface centroids
 *   4. kappa = 0
 *
 * 2D uniform Cartesian (QUAD4) only.
 */
class HeightFunctionCurvature : public GeneralUserObject
{
public:
  static InputParameters validParams();
  HeightFunctionCurvature(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}

protected:
  /// Basilisk height encoding constant
  static constexpr Real HSHIFT = 20.0;
  /// Sentinel for invalid/missing height data
  static constexpr Real HF_NODATA = 1e30;

  /// Decode stored height (strip HSHIFT orientation encoding)
  static Real hfHeight(Real H);
  /// Get orientation from encoded height (0 = full below, 1 = full above)
  static int hfOrientation(Real H);

  /**
   * Compute heights for all cells in a given direction.
   * Runs the two-pass half-column algorithm (j=-1 then j=+1).
   * @param fwd_side QUAD4 side for positive direction (e.g. 2=top for y)
   * @param bck_side QUAD4 side for negative direction (e.g. 0=bottom for y)
   * @param[out] heights Map from element ID to encoded height
   */
  void computeHeightsInDirection(unsigned int fwd_side,
                                 unsigned int bck_side,
                                 std::unordered_map<dof_id_type, Real> & heights);

  /**
   * Half-column computation for one cell (matches Basilisk's half_column).
   * @param elem The cell to compute height for
   * @param j Direction: +1 for forward, -1 for backward
   * @param side QUAD4 side to traverse in direction j
   * @param heights Map storing intermediate/final heights
   */
  void halfColumn(const Elem * elem,
                  int j,
                  unsigned int side,
                  std::unordered_map<dof_id_type, Real> & heights);

  /**
   * Propagate heights along column direction (Basilisk's column_propagation).
   * If a neighbor's height, adjusted by offset, gives a better (closer to 0)
   * estimate, adopt it.
   */
  void columnPropagation(unsigned int fwd_side,
                         unsigned int bck_side,
                         std::unordered_map<dof_id_type, Real> & heights);

  /**
   * Compute curvature from heights in one direction.
   * Uses the 3-point stencil in the differentiation direction.
   * Returns HF_NODATA if the stencil is incomplete or orientations don't match.
   * @param plus_side Side for +1 neighbor in differentiation direction
   * @param minus_side Side for -1 neighbor in differentiation direction
   */
  Real kappaInDirection(const Elem * elem,
                        unsigned int plus_side,
                        unsigned int minus_side,
                        const std::unordered_map<dof_id_type, Real> & heights,
                        Real Delta) const;

  /**
   * Height-function curvature with direction selection and sign flip.
   * Matches Basilisk's height_curvature: sorts directions by |c[+1]-c[-1]|,
   * tries largest first, applies sign correction.
   */
  bool heightCurvature(const Elem * elem, Real Delta, Real & kappa) const;

  /// Check if cell is interfacial (Basilisk convention)
  bool isInterfacial(const Elem * elem) const;

  /**
   * Parabolic fit of PLIC centroids (Basilisk's centroids_curvature_fit).
   */
  bool centroidsCurvatureFit(const Elem * elem,
                             const VectorValue<Real> & n_hat,
                             Real d_global,
                             Real & kappa) const;

  /// Read alpha for a cell
  Real getAlpha(const Elem * elem) const;

  /// The VOF variable
  MooseLinearVariableFV<Real> & _alpha_var;

  /// Linear system for alpha
  LinearSystem & _alpha_system;

  /// PLIC reconstruction (provides interface normals)
  const PLICReconstruction & _plic;

  /// The kappa AuxVariable name
  const VariableName _kappa_name;

  /// Half-extent of height columns (Basilisk default: 4)
  const unsigned int _column_half_extent;

  /// Tolerance for pure/interface cell detection
  const Real _alpha_tol;

  // Cached during execute()
  const NumericVector<Number> * _alpha_solution;
  unsigned int _alpha_sys_num;
  unsigned int _alpha_var_num;

  /// Heights from x-direction integration (columns along x, curvature diff in y)
  std::unordered_map<dof_id_type, Real> _heights_x;
  /// Heights from y-direction integration (columns along y, curvature diff in x)
  std::unordered_map<dof_id_type, Real> _heights_y;
};
