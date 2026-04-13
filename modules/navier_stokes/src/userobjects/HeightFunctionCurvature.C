//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeightFunctionCurvature.h"
#include "PLICReconstruction.h"
#include "MooseLinearVariableFV.h"
#include "FEProblemBase.h"
#include "LinearSystem.h"
#include "MooseMesh.h"
#include "ElemInfo.h"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/elem.h"

registerMooseObject("NavierStokesApp", HeightFunctionCurvature);

// ============================================================================
// ParabolaFit2D — Weighted least-squares parabola in local frame
// Based on Basilisk's parabola.h (Popinet 2009)
// Fits: eta = a[0]*xi^2 + a[1]*xi + a[2]
// where xi is tangent to interface, eta is normal
// ============================================================================

void
ParabolaFit2D::init(const Point & o, const VectorValue<Real> & normal)
{
  origin = o;
  Real mag = normal.norm();
  nx = normal(0) / mag;
  ny = normal(1) / mag;
  tx = ny;
  ty = -nx;
  for (int i = 0; i < 3; i++)
  {
    rhs[i] = 0.0;
    a[i] = 0.0;
    for (int j = 0; j < 3; j++)
      M[i][j] = 0.0;
  }
}

void
ParabolaFit2D::addPoint(const Point & p, Real w)
{
  Real dx = p(0) - origin(0);
  Real dy = p(1) - origin(1);
  Real xi = tx * dx + ty * dy;
  Real eta = nx * dx + ny * dy;

  Real xi2 = xi * xi;
  Real xi3 = xi2 * xi;
  Real xi4 = xi3 * xi;

  M[0][0] += w * xi4;
  M[1][0] += w * xi3;
  M[1][1] += w * xi2;
  M[2][1] += w * xi;
  M[2][2] += w;

  rhs[0] += w * xi2 * eta;
  rhs[1] += w * xi * eta;
  rhs[2] += w * eta;
}

bool
ParabolaFit2D::solve()
{
  M[0][1] = M[1][0];
  M[0][2] = M[1][1];
  M[2][0] = M[1][1];
  M[1][2] = M[2][1];

  Real A[3][3], b[3];
  for (int i = 0; i < 3; i++)
  {
    b[i] = rhs[i];
    for (int j = 0; j < 3; j++)
      A[i][j] = M[i][j];
  }

  for (int col = 0; col < 3; col++)
  {
    int pivot = col;
    Real max_val = std::abs(A[col][col]);
    for (int row = col + 1; row < 3; row++)
    {
      if (std::abs(A[row][col]) > max_val)
      {
        max_val = std::abs(A[row][col]);
        pivot = row;
      }
    }
    if (max_val < 1e-10)
    {
      a[0] = a[1] = a[2] = 0.0;
      return false;
    }
    if (pivot != col)
    {
      std::swap(b[col], b[pivot]);
      for (int j = 0; j < 3; j++)
        std::swap(A[col][j], A[pivot][j]);
    }
    for (int row = col + 1; row < 3; row++)
    {
      Real factor = A[row][col] / A[col][col];
      for (int j = col; j < 3; j++)
        A[row][j] -= factor * A[col][j];
      b[row] -= factor * b[col];
    }
  }

  for (int i = 2; i >= 0; i--)
  {
    a[i] = b[i];
    for (int j = i + 1; j < 3; j++)
      a[i] -= A[i][j] * a[j];
    a[i] /= A[i][i];
  }
  return true;
}

Real
ParabolaFit2D::curvature(Real kappa_max) const
{
  Real denom = std::pow(1.0 + a[1] * a[1], 1.5);
  Real kappa = -2.0 * a[0] / denom;
  return std::max(-kappa_max, std::min(kappa_max, kappa));
}

// ============================================================================
// HeightFunctionCurvature — Basilisk-style height functions
// ============================================================================

InputParameters
HeightFunctionCurvature::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Computes interface curvature using height functions (Popinet 2009). "
      "Heights match Basilisk convention: signed distance from cell center "
      "to interface, normalized by Delta. 2D uniform Cartesian only.");
  params.addRequiredParam<VariableName>("alpha_variable", "The VOF phase fraction variable.");
  params.addRequiredParam<VariableName>("kappa_variable",
                                        "The AuxVariable to write curvature into.");
  params.addRequiredParam<UserObjectName>("plic_reconstruction",
                                          "The PLICReconstruction user object.");
  params.addParam<unsigned int>(
      "column_half_extent", 4, "Half-extent of height columns (Basilisk default: 4).");
  params.addParam<Real>("alpha_tolerance", 1e-6, "Tolerance for pure/interface detection.");
  return params;
}

HeightFunctionCurvature::HeightFunctionCurvature(const InputParameters & params)
  : GeneralUserObject(params),
    _alpha_var(dynamic_cast<MooseLinearVariableFV<Real> &>(
        _fe_problem.getVariable(0, getParam<VariableName>("alpha_variable")))),
    _alpha_system(_fe_problem.getLinearSystem(_alpha_var.sys().number())),
    _plic(getUserObject<PLICReconstruction>("plic_reconstruction")),
    _kappa_name(getParam<VariableName>("kappa_variable")),
    _column_half_extent(getParam<unsigned int>("column_half_extent")),
    _alpha_tol(getParam<Real>("alpha_tolerance")),
    _alpha_solution(nullptr),
    _alpha_sys_num(0),
    _alpha_var_num(0)
{
  _alpha_var.computeCellGradients();
}

Real
HeightFunctionCurvature::hfHeight(Real H)
{
  return H > HSHIFT / 2. ? H - HSHIFT : H < -HSHIFT / 2. ? H + HSHIFT : H;
}

int
HeightFunctionCurvature::hfOrientation(Real H)
{
  return std::abs(H) > HSHIFT / 2. ? 1 : 0;
}

Real
HeightFunctionCurvature::getAlpha(const Elem * elem) const
{
  const auto & ei = _fe_problem.mesh().elemInfo(elem->id());
  const auto alpha_dof = ei.dofIndices()[_alpha_sys_num][_alpha_var_num];
  Real alpha = (*_alpha_solution)(alpha_dof);
  // Clamp to [0,1] with tolerance — Basilisk's VOF is sharp (exact 0/1 in pure cells)
  // but our tanh IC and advection scheme produce diffuse tails. The height-function
  // state machine requires exact 0/1 for transition detection.
  if (alpha < _alpha_tol)
    return 0.0;
  if (alpha > 1.0 - _alpha_tol)
    return 1.0;
  return alpha;
}

// ============================================================================
// Half-column algorithm (matches Basilisk's half_column in heights.h)
//
// For j=-1 (backward pass): computes partial column info, stores intermediate
// state in the heights map (encoded as: 300 = failed, H+100*(1+(S>=1)) = partial,
// or completed height with HSHIFT encoding).
//
// For j=+1 (forward pass): reads intermediate state from j=-1, continues the
// column, combines results keeping the better (smaller |height|) estimate.
// ============================================================================

void
HeightFunctionCurvature::halfColumn(const Elem * elem,
                                    int j,
                                    unsigned int side,
                                    std::unordered_map<dof_id_type, Real> & heights)
{
  const dof_id_type eid = elem->id();

  Real S = getAlpha(elem);
  Real H = S;

  // State from previous (j=-1) pass
  bool state_complete = false;
  Real state_h = 0.0;

  if (j == 1)
  {
    auto it = heights.find(eid);
    if (it != heights.end())
    {
      Real stored = it->second;
      if (stored == 300.0)
      {
        state_complete = true;
        state_h = HF_NODATA;
      }
      else if (stored >= HF_NODATA)
      {
        // Already nodata from a previous direction — skip
        state_complete = true;
        state_h = HF_NODATA;
      }
      else
      {
        int s = static_cast<int>((stored + HSHIFT / 2.0) / 100.0);
        state_h = stored - 100.0 * s;
        int state_s = s - 1;
        if (state_s == -1)
        {
          // Downward pass completed
          state_complete = true;
          // state_h already holds the completed height
        }
        else
        {
          // Continue from partial state
          S = static_cast<Real>(state_s);
          H = state_h;
        }
      }
    }
    else
    {
      // No data from backward pass
      state_complete = true;
      state_h = HF_NODATA;
    }
  }

  // Traverse in direction j
  bool found_complete = false;
  const Elem * curr = elem;
  int i;
  for (i = 1; i <= static_cast<int>(_column_half_extent); ++i)
  {
    const Elem * next = curr->neighbor_ptr(side);
    if (!next)
      break;

    Real ci = getAlpha(next);
    H += ci;

    if (S > 0.0 && S < 1.0)
    {
      S = ci;
      if (ci <= 0.0 || ci >= 1.0)
      {
        H -= i * ci;
        break;
      }
    }
    else if (S >= 1.0 && ci <= 0.0)
    {
      H = (H - 0.5) * j + (j == -1 ? HSHIFT : 0.0);
      found_complete = true;
      break;
    }
    else if (S <= 0.0 && ci >= 1.0)
    {
      H = (i + 0.5 - H) * j + (j == 1 ? HSHIFT : 0.0);
      found_complete = true;
      break;
    }
    else if (S == ci && std::trunc(H) != H)
      break;

    curr = next;
  }

  if (j == -1)
  {
    Real alpha_cell = getAlpha(elem);
    if (!found_complete &&
        ((alpha_cell <= 0.0 || alpha_cell >= 1.0) || (S > 0.0 && S < 1.0)))
      heights[eid] = 300.0;
    else if (found_complete)
      heights[eid] = H;
    else
      heights[eid] = H + 100.0 * (1.0 + (S >= 1.0 ? 1.0 : 0.0));
  }
  else // j == +1
  {
    if (!state_complete ||
        (found_complete && std::abs(hfHeight(H)) < std::abs(hfHeight(state_h))))
    {
      state_complete = found_complete;
      state_h = H;
    }

    if (!state_complete)
      heights[eid] = HF_NODATA;
    else
      heights[eid] = (state_h > 1e10 ? HF_NODATA : state_h);
  }
}

void
HeightFunctionCurvature::computeHeightsInDirection(
    unsigned int fwd_side,
    unsigned int bck_side,
    std::unordered_map<dof_id_type, Real> & heights)
{
  auto & mesh = _fe_problem.mesh();
  heights.clear();

  // Pass 1: backward (j = -1)
  for (const auto * elem_info_ptr : mesh.elemInfoVector())
    halfColumn(elem_info_ptr->elem(), -1, bck_side, heights);

  // Pass 2: forward (j = +1)
  for (const auto * elem_info_ptr : mesh.elemInfoVector())
    halfColumn(elem_info_ptr->elem(), +1, fwd_side, heights);
}

void
HeightFunctionCurvature::columnPropagation(
    unsigned int fwd_side,
    unsigned int bck_side,
    std::unordered_map<dof_id_type, Real> & heights)
{
  auto & mesh = _fe_problem.mesh();

  for (const auto * elem_info_ptr : mesh.elemInfoVector())
  {
    const auto * elem = elem_info_ptr->elem();
    const dof_id_type eid = elem->id();

    auto it = heights.find(eid);
    Real h_curr = (it != heights.end()) ? it->second : HF_NODATA;

    // Check neighbors at offsets -2 to +2 along the column direction
    for (int offset = -2; offset <= 2; ++offset)
    {
      if (offset == 0)
        continue;

      // Navigate to neighbor at |offset| steps
      const Elem * nb = elem;
      unsigned int nav_side = (offset > 0) ? fwd_side : bck_side;
      int abs_off = std::abs(offset);
      bool valid = true;
      for (int k = 0; k < abs_off; ++k)
      {
        if (!nb)
        {
          valid = false;
          break;
        }
        nb = nb->neighbor_ptr(nav_side);
      }
      if (!valid || !nb)
        continue;

      auto nb_it = heights.find(nb->id());
      if (nb_it == heights.end())
        continue;

      Real h_nb = nb_it->second;
      if (h_nb >= HF_NODATA)
        continue;

      Real height_nb = hfHeight(h_nb);
      if (std::abs(height_nb) <= 3.5 &&
          std::abs(height_nb + offset) < std::abs(hfHeight(h_curr)))
      {
        h_curr = h_nb + offset;
        heights[eid] = h_curr;
      }
    }
  }
}

Real
HeightFunctionCurvature::kappaInDirection(
    const Elem * elem,
    unsigned int plus_side,
    unsigned int minus_side,
    const std::unordered_map<dof_id_type, Real> & heights,
    Real Delta) const
{
  // Center height
  auto it_c = heights.find(elem->id());
  if (it_c == heights.end() || it_c->second >= HF_NODATA)
    return HF_NODATA;

  Real h_c = it_c->second;
  int ori = hfOrientation(h_c);

  // +1 neighbor in differentiation direction
  const Elem * plus = elem->neighbor_ptr(plus_side);
  if (!plus)
    return HF_NODATA;
  auto it_p = heights.find(plus->id());
  if (it_p == heights.end() || it_p->second >= HF_NODATA ||
      hfOrientation(it_p->second) != ori)
    return HF_NODATA;

  // -1 neighbor in differentiation direction
  const Elem * minus = elem->neighbor_ptr(minus_side);
  if (!minus)
    return HF_NODATA;
  auto it_m = heights.find(minus->id());
  if (it_m == heights.end() || it_m->second >= HF_NODATA ||
      hfOrientation(it_m->second) != ori)
    return HF_NODATA;

  Real hp = hfHeight(it_p->second);
  Real hm = hfHeight(it_m->second);
  Real h0 = hfHeight(h_c);

  // Basilisk's curvature formula (heights.h / curvature.h):
  // hx = (h[+1] - h[-1]) / 2        (dimensionless, grid-index spacing = 1)
  // hxx = (h[+1] + h[-1] - 2*h[0]) / Delta  (physical second derivative)
  // kappa = hxx / (1 + hx^2)^{3/2}
  Real hx = (hp - hm) / 2.0;
  Real hxx = (hp + hm - 2.0 * h0) / Delta;
  Real denom = std::pow(1.0 + hx * hx, 1.5);

  return hxx / denom;
}

bool
HeightFunctionCurvature::heightCurvature(const Elem * elem,
                                         Real Delta,
                                         Real & kappa) const
{
  // Normal estimation: n = c[+1] - c[-1] per direction (Basilisk convention)
  // QUAD4 sides: 0=bottom, 1=right, 2=top, 3=left
  const Elem * right = elem->neighbor_ptr(1);
  const Elem * left = elem->neighbor_ptr(3);
  const Elem * top = elem->neighbor_ptr(2);
  const Elem * bottom = elem->neighbor_ptr(0);

  Real n_x = (right ? getAlpha(right) : getAlpha(elem)) -
             (left ? getAlpha(left) : getAlpha(elem));
  Real n_y = (top ? getAlpha(top) : getAlpha(elem)) -
             (bottom ? getAlpha(bottom) : getAlpha(elem));

  // Direction info: {normal component, height map, diff +side, diff -side}
  // kappa_x: from _heights_x (columns along x), differentiate in y (top=2, bottom=0)
  // kappa_y: from _heights_y (columns along y), differentiate in x (right=1, left=3)
  struct DirInfo
  {
    Real n;
    const std::unordered_map<dof_id_type, Real> * heights;
    unsigned int plus_side;
    unsigned int minus_side;
  };

  DirInfo dirs[2] = {{n_x, &_heights_x, 2, 0},   // x-heights, diff in y
                     {n_y, &_heights_y, 1, 3}};   // y-heights, diff in x

  // Sort: largest |n| first (Basilisk sorts dimensions by normal magnitude)
  if (std::abs(dirs[0].n) < std::abs(dirs[1].n))
    std::swap(dirs[0], dirs[1]);

  for (int d = 0; d < 2; ++d)
  {
    Real k = kappaInDirection(
        elem, dirs[d].plus_side, dirs[d].minus_side, *dirs[d].heights, Delta);
    if (k < HF_NODATA)
    {
      // Sign flip based on normal direction (Basilisk: if n < 0, kappa = -kappa)
      if (dirs[d].n < 0.0)
        k = -k;

      // Clamp: |kappa| <= 1/Delta
      if (std::abs(k) > 1.0 / Delta)
        k = (k > 0 ? 1.0 : -1.0) / Delta;

      kappa = k;
      return true;
    }
  }

  return false;
}

bool
HeightFunctionCurvature::isInterfacial(const Elem * elem) const
{
  Real c = getAlpha(elem);
  if (c > 0.0 && c < 1.0)
    return true;

  for (unsigned int s = 0; s < elem->n_sides(); ++s)
  {
    const Elem * nb = elem->neighbor_ptr(s);
    if (!nb)
      continue;
    Real cn = getAlpha(nb);
    if (c >= 1.0 && cn <= 0.0)
      return true;
    if (c <= 0.0 && cn >= 1.0)
      return true;
  }
  return false;
}

bool
HeightFunctionCurvature::centroidsCurvatureFit(const Elem * elem,
                                               const VectorValue<Real> & n_hat,
                                               Real d_global,
                                               Real & kappa) const
{
  const Point cell_center = elem->vertex_average();
  const Real dist = n_hat(0) * cell_center(0) + n_hat(1) * cell_center(1) - d_global;
  Point centroid(cell_center(0) - dist * n_hat(0),
                 cell_center(1) - dist * n_hat(1),
                 0.0);

  const Real dx = std::abs(elem->point(2)(0) - elem->point(0)(0));

  ParabolaFit2D fit;
  fit.init(centroid, n_hat);
  fit.addPoint(centroid, 0.1);

  int n_points = 1;
  for (unsigned int s = 0; s < elem->n_sides(); ++s)
  {
    const Elem * nb = elem->neighbor_ptr(s);
    if (!nb)
      continue;

    if (_plic.hasPlane(nb->id()))
    {
      const auto & nb_plane = _plic.getPlane(nb->id());
      const Point nb_center = nb->vertex_average();
      const Real nb_dist =
          nb_plane.n_hat(0) * nb_center(0) + nb_plane.n_hat(1) * nb_center(1) - nb_plane.d;
      Point nb_centroid(nb_center(0) - nb_dist * nb_plane.n_hat(0),
                        nb_center(1) - nb_dist * nb_plane.n_hat(1),
                        0.0);
      fit.addPoint(nb_centroid, 1.0);
      n_points++;
    }

    // Diagonal neighbors
    for (unsigned int s2 = 0; s2 < nb->n_sides(); ++s2)
    {
      const Elem * nb2 = nb->neighbor_ptr(s2);
      if (!nb2 || nb2 == elem)
        continue;
      bool is_face_neighbor = false;
      for (unsigned int sf = 0; sf < elem->n_sides(); ++sf)
        if (elem->neighbor_ptr(sf) == nb2)
        {
          is_face_neighbor = true;
          break;
        }
      if (is_face_neighbor)
        continue;

      if (_plic.hasPlane(nb2->id()))
      {
        const auto & nb2_plane = _plic.getPlane(nb2->id());
        const Point nb2_center = nb2->vertex_average();
        const Real nb2_dist = nb2_plane.n_hat(0) * nb2_center(0) +
                              nb2_plane.n_hat(1) * nb2_center(1) - nb2_plane.d;
        Point nb2_centroid(nb2_center(0) - nb2_dist * nb2_plane.n_hat(0),
                           nb2_center(1) - nb2_dist * nb2_plane.n_hat(1),
                           0.0);
        fit.addPoint(nb2_centroid, 1.0);
        n_points++;
      }
    }
  }

  if (n_points < 3)
    return false;

  if (!fit.solve())
    return false;

  kappa = fit.curvature(1.0 / dx);
  return true;
}

void
HeightFunctionCurvature::execute()
{
  auto & mesh = _fe_problem.mesh();

  // Cache alpha solution access
  const auto & alpha_sys =
      libMesh::cast_ref<libMesh::LinearImplicitSystem &>(_alpha_var.sys().system());
  _alpha_solution = alpha_sys.solution.get();
  _alpha_sys_num = _alpha_var.sys().number();
  _alpha_var_num = _alpha_var.number();

  // Ensure alpha gradients are up to date
  _alpha_system.computeGradients();

  // ---- Step 1: Compute heights for all cells ----
  // x-direction: columns along x (fwd=right=1, bck=left=3)
  computeHeightsInDirection(1, 3, _heights_x);
  columnPropagation(1, 3, _heights_x);

  // y-direction: columns along y (fwd=top=2, bck=bottom=0)
  computeHeightsInDirection(2, 0, _heights_y);
  columnPropagation(2, 0, _heights_y);

  // ---- Step 2: Compute curvature ----
  auto & kappa_var =
      dynamic_cast<MooseVariableFVReal &>(_fe_problem.getVariable(0, _kappa_name));
  auto & kappa_sys = kappa_var.sys().system();
  auto & kappa_solution = *kappa_sys.solution;
  auto & petsc_kappa = dynamic_cast<PetscVector<Number> &>(kappa_solution);
  PetscScalar * kappa_array = petsc_kappa.get_array();

  const auto kappa_sys_num = kappa_var.sys().number();
  const auto kappa_var_num = kappa_var.number();

  // Get cell size (uniform mesh)
  const auto * first_elem = mesh.elemInfoVector()[0]->elem();
  const Real Delta = std::abs(first_elem->point(2)(0) - first_elem->point(0)(0));

  constexpr Real NEEDS_FALLBACK = std::numeric_limits<Real>::max();

  unsigned int n_interface = 0, n_hf = 0;
  Real kappa_min = 1e30, kappa_max = -1e30;

  // ---- Pass 1: Height-function curvature ----
  for (const auto * elem_info_ptr : mesh.elemInfoVector())
  {
    const auto & ei = *elem_info_ptr;
    const auto kappa_dof = ei.dofIndices()[kappa_sys_num][kappa_var_num];

    if (kappa_dof < petsc_kappa.first_local_index() ||
        kappa_dof >= petsc_kappa.last_local_index())
      continue;

    const auto local_id = kappa_dof - petsc_kappa.first_local_index();

    if (!isInterfacial(ei.elem()))
    {
      kappa_array[local_id] = 0.0;
      continue;
    }

    n_interface++;
    Real kappa_val = 0.0;

    if (heightCurvature(ei.elem(), Delta, kappa_val))
    {
      // Negate: HF produces Basilisk convention (kappa = +1/R for convex),
      // but MOOSE's CSF kernel expects div(n_hat) convention (kappa = -1/R).
      kappa_array[local_id] = -kappa_val;
      n_hf++;
      if (kappa_val < kappa_min) kappa_min = kappa_val;
      if (kappa_val > kappa_max) kappa_max = kappa_val;
    }
    else
      kappa_array[local_id] = NEEDS_FALLBACK;
  }

  // ---- Pass 2: Neighbor averaging ----
  for (const auto * elem_info_ptr : mesh.elemInfoVector())
  {
    const auto & ei = *elem_info_ptr;
    const auto kappa_dof = ei.dofIndices()[kappa_sys_num][kappa_var_num];

    if (kappa_dof < petsc_kappa.first_local_index() ||
        kappa_dof >= petsc_kappa.last_local_index())
      continue;

    const auto local_id = kappa_dof - petsc_kappa.first_local_index();

    if (kappa_array[local_id] != NEEDS_FALLBACK)
      continue;

    Real sum = 0.0;
    unsigned int count = 0;
    const auto * elem = ei.elem();

    // Face and diagonal neighbors
    for (unsigned int s = 0; s < elem->n_sides(); ++s)
    {
      const Elem * nb = elem->neighbor_ptr(s);
      if (!nb)
        continue;
      const auto & nb_ei = mesh.elemInfo(nb->id());
      const auto nb_dof = nb_ei.dofIndices()[kappa_sys_num][kappa_var_num];
      if (nb_dof < petsc_kappa.first_local_index() ||
          nb_dof >= petsc_kappa.last_local_index())
        continue;
      const auto nb_local = nb_dof - petsc_kappa.first_local_index();
      const Real nb_kappa = kappa_array[nb_local];
      if (nb_kappa != NEEDS_FALLBACK && nb_kappa != 0.0)
      {
        sum += nb_kappa;
        count++;
      }

      for (unsigned int s2 = 0; s2 < nb->n_sides(); ++s2)
      {
        const Elem * nb2 = nb->neighbor_ptr(s2);
        if (!nb2 || nb2 == elem)
          continue;
        bool is_face = false;
        for (unsigned int sf = 0; sf < elem->n_sides(); ++sf)
          if (elem->neighbor_ptr(sf) == nb2)
          {
            is_face = true;
            break;
          }
        if (is_face)
          continue;
        const auto & nb2_ei = mesh.elemInfo(nb2->id());
        const auto nb2_dof = nb2_ei.dofIndices()[kappa_sys_num][kappa_var_num];
        if (nb2_dof < petsc_kappa.first_local_index() ||
            nb2_dof >= petsc_kappa.last_local_index())
          continue;
        const auto nb2_local = nb2_dof - petsc_kappa.first_local_index();
        const Real nb2_kappa = kappa_array[nb2_local];
        if (nb2_kappa != NEEDS_FALLBACK && nb2_kappa != 0.0)
        {
          sum += nb2_kappa;
          count++;
        }
      }
    }

    if (count > 0)
    {
      kappa_array[local_id] = sum / count;
    }
  }

  // ---- Pass 3: Centroid parabolic fit ----
  for (const auto * elem_info_ptr : mesh.elemInfoVector())
  {
    const auto & ei = *elem_info_ptr;
    const auto kappa_dof = ei.dofIndices()[kappa_sys_num][kappa_var_num];

    if (kappa_dof < petsc_kappa.first_local_index() ||
        kappa_dof >= petsc_kappa.last_local_index())
      continue;

    const auto local_id = kappa_dof - petsc_kappa.first_local_index();

    if (kappa_array[local_id] != NEEDS_FALLBACK)
      continue;

    const auto elem_id = ei.elem()->id();
    if (!_plic.hasPlane(elem_id))
    {
      kappa_array[local_id] = 0.0;
      continue;
    }

    const auto & plane = _plic.getPlane(elem_id);
    Real kappa_val = 0.0;

    if (centroidsCurvatureFit(ei.elem(), plane.n_hat, plane.d, kappa_val))
      kappa_array[local_id] = -kappa_val;
    else
      kappa_array[local_id] = 0.0;
  }

  mooseWarning("HF curvature: ", n_interface, " interface, ", n_hf, "/", n_interface,
               " HF, kappa=[", kappa_min, ",", kappa_max, "] Delta=", Delta);

  petsc_kappa.restore_array();
  kappa_solution.close();

  // Sync current_local_solution so functor evaluations (e.g. BalancedForceSurfaceTension)
  // see the updated kappa values
  kappa_sys.update();
}
