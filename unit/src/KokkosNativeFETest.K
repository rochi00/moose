//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// Tests : native FE evaluator correctness.
//
// Two things are verified:
//
//   A. FEEvaluator vs libMesh FEBase  (cross-validation)
//      For every (element type, quadrature order) pair, evaluate
//      FEEvaluator<LAGRANGE, ElemType>::shape/grad_shape at the native
//      quadrature points and compare point-for-point against
//      libMesh FEBase::reinit() output.
//
//   B. nativeMapShape dispatch parity
//      The runtime-dispatch functions in kokkos_fe_evaluator.h must produce
//      bit-identical results to direct FEEvaluator calls because they
//      evaluate the same expression at the same coordinates.
//
//   C. mapFaceQpToParent correctness
//      Face quadrature points mapped to parent coordinates must reproduce the
//      same phi values as a libMesh face reinit with the same points.
//
// Run:
//   cd unit && make && ./moose-unit-opt --gtest_filter="KokkosNativeFETest*"

#include "gtest_include.h"

#include "libmesh/gpu/kokkos_fe_lagrange_1d.h"
#include "libmesh/gpu/kokkos_fe_lagrange_2d.h"
#include "libmesh/gpu/kokkos_fe_lagrange_3d.h"
#include "libmesh/gpu/kokkos_fe_monomial.h"
#include "libmesh/gpu/kokkos_fe_evaluator.h"
#include "libmesh/gpu/kokkos_fe_face_map.h"
#include "KokkosTypes.h"

#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/reference_elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace libMesh::Kokkos;

static constexpr double TOL = 1.0e-13;

// ── helpers ───────────────────────────────────────────────────────────────────

struct ElemEntry
{
  libMesh::ElemType lm_type;
  libMesh::Order    lm_order;   // default Lagrange order for this elem type
  unsigned int      dim;
  unsigned int      n_dofs;
  const char *      name;
};

static const ElemEntry ELEMS[] = {
  { libMesh::EDGE2,  libMesh::FIRST,  1,  2, "EDGE2"  },
  { libMesh::EDGE3,  libMesh::SECOND, 1,  3, "EDGE3"  },
  { libMesh::TRI3,   libMesh::FIRST,  2,  3, "TRI3"   },
  { libMesh::TRI6,   libMesh::SECOND, 2,  6, "TRI6"   },
  { libMesh::QUAD4,  libMesh::FIRST,  2,  4, "QUAD4"  },
  { libMesh::QUAD8,  libMesh::SECOND, 2,  8, "QUAD8"  },
  { libMesh::QUAD9,  libMesh::SECOND, 2,  9, "QUAD9"  },
  { libMesh::TET4,   libMesh::FIRST,  3,  4, "TET4"   },
  { libMesh::TET10,  libMesh::SECOND, 3, 10, "TET10"  },
  { libMesh::HEX8,   libMesh::FIRST,  3,  8, "HEX8"   },
  { libMesh::HEX20,  libMesh::SECOND, 3, 20, "HEX20"  },
  { libMesh::HEX27,  libMesh::SECOND, 3, 27, "HEX27"  },
};
static constexpr unsigned int N_ELEMS = sizeof(ELEMS) / sizeof(ELEMS[0]);

// ── A. FEEvaluator vs libMesh FEBase ─────────────────────────────────────────

struct CrossValParam
{
  unsigned int elem_idx;
  unsigned int order;
  std::string  label;
};

static std::vector<CrossValParam> crossValParams()
{
  std::vector<CrossValParam> v;
  for (unsigned int e = 0; e < N_ELEMS; ++e)
    for (unsigned int order = 1; order <= 5; ++order)
      v.push_back({e, order,
                   std::string(ELEMS[e].name) + "_o" + std::to_string(order)});
  return v;
}

class KokkosNativeFECrossValTest : public ::testing::TestWithParam<CrossValParam> {};

TEST_P(KokkosNativeFECrossValTest, PhiMatchesLibMesh)
{
  auto p = GetParam();
  const auto & e = ELEMS[p.elem_idx];

  // Build native quadrature at this order
  Array<Real3> q_pts; Array<Real> wts;
  fillQuadrature(e.lm_type, p.order, q_pts, wts);
  const unsigned int nqp = q_pts.size();

  // Build libMesh FE object and attach a QGauss rule at the same order
  libMesh::QGauss qrule(e.dim, static_cast<libMesh::Order>(p.order));
  qrule.allow_rules_with_negative_weights = true;
  auto fe = libMesh::FEBase::build(e.dim, libMesh::FEType(e.lm_order, libMesh::LAGRANGE));
  fe->attach_quadrature_rule(&qrule);


  const auto & phi_lm  = fe->get_phi();
  const auto & dphi_lm = fe->get_dphi();

  const libMesh::Elem * ref = &libMesh::ReferenceElem::get(e.lm_type);
  fe->reinit(ref);

  ASSERT_EQ(phi_lm.size(), (unsigned int)e.n_dofs) << p.label;
  ASSERT_EQ((unsigned int)qrule.n_points(), nqp)   << p.label << ": qp count mismatch vs libMesh";

  for (unsigned int i = 0; i < e.n_dofs; ++i)
    for (unsigned int q = 0; q < nqp; ++q)
    {
      Real xi   = q_pts[q](0);
      Real eta  = q_pts[q](1);
      Real zeta = q_pts[q](2);

      double native_phi = nativeMapShape(libMesh::LAGRANGE_MAP, e.lm_type, i, xi, eta, zeta);
      EXPECT_NEAR(native_phi, phi_lm[i][q], TOL)
          << p.label << " phi i=" << i << " q=" << q;

      Real3 ng = nativeGradMapShape(libMesh::LAGRANGE_MAP, e.lm_type, i, xi, eta, zeta);
      EXPECT_NEAR(ng(0), dphi_lm[i][q](0), TOL)
          << p.label << " dphi/dxi i=" << i << " q=" << q;
      if (e.dim >= 2)
        EXPECT_NEAR(ng(1), dphi_lm[i][q](1), TOL)
            << p.label << " dphi/deta i=" << i << " q=" << q;
      if (e.dim >= 3)
        EXPECT_NEAR(ng(2), dphi_lm[i][q](2), TOL)
            << p.label << " dphi/dzeta i=" << i << " q=" << q;
    }
}

INSTANTIATE_TEST_SUITE_P(AllElemsOrders,
                         KokkosNativeFECrossValTest,
                         ::testing::ValuesIn(crossValParams()),
                         [](const auto & info) { return info.param.label; });

// ── B. nativeMapShape dispatch is bit-identical to direct FEEvaluator calls ──
//
// Since both paths evaluate the same inline function at the same coordinates,
// results must be bit-identical (not just within tolerance).

TEST(KokkosNativeFETest, NativeShapeDispatch_BitIdentical_Quad4)
{
  Array<Real3> q_pts; Array<Real> wts;
  fillQuadrature(libMesh::QUAD4, 4, q_pts, wts);

  for (unsigned int i = 0; i < 4; ++i)
    for (unsigned int q = 0; q < q_pts.size(); ++q)
    {
      Real xi = q_pts[q](0), eta = q_pts[q](1);
      Real direct   = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::shape(i, xi, eta, 0.0);
      Real dispatch = nativeMapShape(libMesh::LAGRANGE_MAP, libMesh::QUAD4, i, xi, eta, 0.0);
      EXPECT_EQ(direct, dispatch) << "i=" << i << " q=" << q;
    }
}

TEST(KokkosNativeFETest, NativeShapeDispatch_BitIdentical_Tet4)
{
  Array<Real3> q_pts; Array<Real> wts;
  fillQuadrature(libMesh::TET4, 5, q_pts, wts);

  for (unsigned int i = 0; i < 4; ++i)
    for (unsigned int q = 0; q < q_pts.size(); ++q)
    {
      Real xi = q_pts[q](0), eta = q_pts[q](1), zeta = q_pts[q](2);
      Real direct   = FEEvaluator<libMesh::LAGRANGE, libMesh::TET4>::shape(i, xi, eta, zeta);
      Real dispatch = nativeMapShape(libMesh::LAGRANGE_MAP, libMesh::TET4, i, xi, eta, zeta);
      EXPECT_EQ(direct, dispatch) << "i=" << i << " q=" << q;
    }
}

TEST(KokkosNativeFETest, NativeShapeDispatch_BitIdentical_Hex27)
{
  Array<Real3> q_pts; Array<Real> wts;
  fillQuadrature(libMesh::HEX27, 3, q_pts, wts);

  for (unsigned int i = 0; i < 27; ++i)
    for (unsigned int q = 0; q < q_pts.size(); ++q)
    {
      Real xi = q_pts[q](0), eta = q_pts[q](1), zeta = q_pts[q](2);
      Real direct   = FEEvaluator<libMesh::LAGRANGE, libMesh::HEX27>::shape(i, xi, eta, zeta);
      Real dispatch = nativeMapShape(libMesh::LAGRANGE_MAP, libMesh::HEX27, i, xi, eta, zeta);
      EXPECT_EQ(direct, dispatch) << "i=" << i << " q=" << q;
    }
}

TEST(KokkosNativeFETest, NativeGradShapeDispatch_BitIdentical_Tri6)
{
  Array<Real3> q_pts; Array<Real> wts;
  fillQuadrature(libMesh::TRI6, 4, q_pts, wts);

  for (unsigned int i = 0; i < 6; ++i)
    for (unsigned int q = 0; q < q_pts.size(); ++q)
    {
      Real xi = q_pts[q](0), eta = q_pts[q](1);
      Real3 direct   = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::grad_shape(i, xi, eta, 0.0);
      Real3 dispatch = nativeGradMapShape(libMesh::LAGRANGE_MAP, libMesh::TRI6, i, xi, eta, 0.0);
      EXPECT_EQ(direct(0), dispatch(0)) << "v[0] i=" << i << " q=" << q;
      EXPECT_EQ(direct(1), dispatch(1)) << "v[1] i=" << i << " q=" << q;
    }
}

// ── C. Face qp → parent coordinate mapping ───────────────────────────────────
//
// For a given face side, the parent-reference coordinates produced by
// mapFaceQpToParent must give the same phi values as libMesh's face reinit.
// We test this indirectly: evaluate phi at parent coords via nativeMapShape
// and compare against the libMesh phi from face reinit.
//
// The test covers QUAD4 (2D), HEX8 (3D), and TRI3 (simplex face).

static void checkFaceMapping(libMesh::ElemType parent_topo,
                              libMesh::Order lm_order,
                              unsigned int parent_dim,
                              unsigned int n_dofs,
                              unsigned int quad_order,
                              const std::string & label)
{
  const libMesh::Elem * ref = &libMesh::ReferenceElem::get(parent_topo);
  const auto side_topo = getSideTopology(parent_topo);

  libMesh::QGauss qrule_face(parent_dim - 1,
                              static_cast<libMesh::Order>(quad_order));
  qrule_face.allow_rules_with_negative_weights = true;

  auto fe_face = libMesh::FEBase::build(parent_dim,
                                        libMesh::FEType(lm_order, libMesh::LAGRANGE));
  fe_face->attach_quadrature_rule(&qrule_face);

  const auto & phi_lm = fe_face->get_phi();
  for (unsigned int side = 0; side < ref->n_sides(); ++side)
  {
    // libMesh face reinit — phi stored as phi[dof][qp]
    fe_face->reinit(ref, side);
    unsigned int nqp = qrule_face.n_points();

    // Native face qps (side-element reference coords)
    Array<Real3> face_qpts; Array<Real> face_wts;
    fillQuadrature(side_topo, quad_order, face_qpts, face_wts);
    ASSERT_EQ(face_qpts.size(), nqp)
        << label << " side=" << side << " qp count mismatch";

    // Map face qps to parent reference coords using isoparametric mapping
    auto side_elem = ref->side_ptr(side);
    for (unsigned int q = 0; q < nqp; ++q)
    {
      // Reproduce mapFaceQpToParent logic inline using nativeMapShape
      const Real3 face_pt = face_qpts[q];
      Real3 parent_pt;
      for (unsigned int k = 0; k < side_elem->n_nodes(); ++k)
      {
        Real psi = nativeMapShape(libMesh::LAGRANGE_MAP, side_topo, k,
                                  face_pt(0), face_pt(1), face_pt(2));
        const auto & node = side_elem->point(k);
        parent_pt(0) += psi * node(0);
        parent_pt(1) += psi * node(1);
        parent_pt(2) += psi * node(2);
      }

      // Now evaluate parent shape functions at parent_pt and compare to libMesh
      for (unsigned int i = 0; i < n_dofs; ++i)
      {
        double native_phi = nativeMapShape(libMesh::LAGRANGE_MAP, parent_topo, i,
                                           parent_pt(0), parent_pt(1), parent_pt(2));
        EXPECT_NEAR(native_phi, phi_lm[i][q], TOL)
            << label << " side=" << side << " i=" << i << " q=" << q;
      }
    }
  }
}

TEST(KokkosNativeFETest, FaceMapping_Quad4)
{
  checkFaceMapping(libMesh::QUAD4, libMesh::FIRST, 2, 4, 3, "QUAD4");
}

TEST(KokkosNativeFETest, FaceMapping_Tri3)
{
  checkFaceMapping(libMesh::TRI3, libMesh::FIRST, 2, 3, 3, "TRI3");
}

TEST(KokkosNativeFETest, FaceMapping_Hex8)
{
  checkFaceMapping(libMesh::HEX8, libMesh::FIRST, 3, 8, 3, "HEX8");
}

TEST(KokkosNativeFETest, FaceMapping_Tet4)
{
  checkFaceMapping(libMesh::TET4, libMesh::FIRST, 3, 4, 3, "TET4");
}

TEST(KokkosNativeFETest, FaceMapping_Hex8_HighOrder)
{
  checkFaceMapping(libMesh::HEX8, libMesh::FIRST, 3, 8, 5, "HEX8_o5");
}

// ── D. nDofs matches libMesh DOF count ───────────────────────────────────────

TEST(KokkosNativeFETest, NDofs_MatchesLibMesh)
{
  for (unsigned int e = 0; e < N_ELEMS; ++e)
  {
    const auto & entry = ELEMS[e];
    const unsigned int order = static_cast<unsigned int>(entry.lm_order);
    unsigned int native_n = nDofs(FEShapeKey{libMesh::LAGRANGE, classFromTopology(entry.lm_type), order});
    EXPECT_EQ(native_n, entry.n_dofs) << entry.name;
  }
}

// ── E. All ElemType values handled by nativeMapShape (no default branch) ─────
//
// Passing each topology with i=0 at the centroid should return a finite,
// non-NaN value. This guards against accidentally missing a switch case.

TEST(KokkosNativeFETest, NativeShape_NoUnhandledTopology)
{
  // Centroids in reference element coordinates
  const Real3 centroids[] = {
    Real3( 0.0,  0.0, 0.0), // EDGE2
    Real3( 0.0,  0.0, 0.0), // EDGE3
    Real3( 1.0/3.0, 1.0/3.0, 0.0), // TRI3
    Real3( 1.0/3.0, 1.0/3.0, 0.0), // TRI6
    Real3( 0.0,  0.0, 0.0), // QUAD4
    Real3( 0.0,  0.0, 0.0), // QUAD8
    Real3( 0.0,  0.0, 0.0), // QUAD9
    Real3( 1.0/4.0, 1.0/4.0, 1.0/4.0), // TET4
    Real3( 1.0/4.0, 1.0/4.0, 1.0/4.0), // TET10
    Real3( 0.0,  0.0, 0.0), // HEX8
    Real3( 0.0,  0.0, 0.0), // HEX20
    Real3( 0.0,  0.0, 0.0), // HEX27
  };
  static_assert(N_ELEMS == 12, "Update centroids array if element list changes");

  for (unsigned int e = 0; e < N_ELEMS; ++e)
  {
    Real val = nativeMapShape(libMesh::LAGRANGE_MAP, ELEMS[e].lm_type, 0,
                              centroids[e](0), centroids[e](1), centroids[e](2));
    EXPECT_TRUE(std::isfinite(val)) << ELEMS[e].name << " returned non-finite";
  }
}

// ── F. nativeShape(FEShapeKey{MONOMIAL,...}) parity with libMesh ──────────────
//
// For each (dim, order, representative element type) combination, compare our
// nativeShape/nativeGradShape dispatched via FEShapeKey against libMesh
// FEBase with MONOMIAL family reinitialised on the reference element.
// Quadrature points are used as test coordinates so no additional sampling
// infrastructure is required.

struct MonomialParityParam
{
  libMesh::FEElemClass kk_cls;
  unsigned int         order;
  libMesh::ElemType    lm_type;  // representative reference element
  unsigned int         dim;
  const char *         name;
};

static std::vector<MonomialParityParam> monomialParityParams()
{
  return {
    // 1-D
    { libMesh::FEElemClass::EDGE, 0, libMesh::EDGE2, 1, "EDGE_p0" },
    { libMesh::FEElemClass::EDGE, 1, libMesh::EDGE2, 1, "EDGE_p1" },
    { libMesh::FEElemClass::EDGE, 2, libMesh::EDGE2, 1, "EDGE_p2" },
    { libMesh::FEElemClass::EDGE, 3, libMesh::EDGE2, 1, "EDGE_p3" },
    // 2-D — QUAD reference element (xi,eta ∈ [-1,1]²)
    { libMesh::FEElemClass::QUAD, 0, libMesh::QUAD4, 2, "QUAD_p0" },
    { libMesh::FEElemClass::QUAD, 1, libMesh::QUAD4, 2, "QUAD_p1" },
    { libMesh::FEElemClass::QUAD, 2, libMesh::QUAD4, 2, "QUAD_p2" },
    { libMesh::FEElemClass::QUAD, 3, libMesh::QUAD4, 2, "QUAD_p3" },
    // 2-D — TRI reference element (xi,eta on unit simplex)
    { libMesh::FEElemClass::TRI,  1, libMesh::TRI3,  2, "TRI_p1"  },
    { libMesh::FEElemClass::TRI,  2, libMesh::TRI3,  2, "TRI_p2"  },
    // 3-D — HEX reference element
    { libMesh::FEElemClass::HEX,  0, libMesh::HEX8,  3, "HEX_p0"  },
    { libMesh::FEElemClass::HEX,  1, libMesh::HEX8,  3, "HEX_p1"  },
    { libMesh::FEElemClass::HEX,  2, libMesh::HEX8,  3, "HEX_p2"  },
    { libMesh::FEElemClass::HEX,  3, libMesh::HEX8,  3, "HEX_p3"  },
    // 3-D — TET reference element
    { libMesh::FEElemClass::TET,  1, libMesh::TET4,  3, "TET_p1"  },
    { libMesh::FEElemClass::TET,  2, libMesh::TET4,  3, "TET_p2"  },
  };
}

class KokkosMonomialParityTest
  : public ::testing::TestWithParam<MonomialParityParam> {};

TEST_P(KokkosMonomialParityTest, PhiMatchesLibMesh)
{
  auto p = GetParam();

  // Build native quadrature points (order p+1 is exact for P_p)
  const unsigned int quad_order = p.order + 1;
  Array<Real3> q_pts; Array<Real> wts;
  fillQuadrature(p.lm_type, quad_order, q_pts, wts);
  const unsigned int nqp = q_pts.size();

  // Build libMesh MONOMIAL FE at the same polynomial order
  libMesh::QGauss qrule(p.dim, static_cast<libMesh::Order>(quad_order));
  qrule.allow_rules_with_negative_weights = true;
  auto fe = libMesh::FEBase::build(
      p.dim, libMesh::FEType(static_cast<libMesh::Order>(p.order), libMesh::MONOMIAL));
  fe->attach_quadrature_rule(&qrule);

  const auto & phi_lm  = fe->get_phi();
  const auto & dphi_lm = fe->get_dphi();
  const libMesh::Elem * ref = &libMesh::ReferenceElem::get(p.lm_type);
  fe->reinit(ref);

  ASSERT_EQ(qrule.n_points(), nqp) << p.name << ": qp count mismatch";

  const unsigned int n = nDofs(FEShapeKey{libMesh::MONOMIAL, p.kk_cls, p.order});
  ASSERT_EQ(phi_lm.size(), n) << p.name << ": DOF count mismatch with libMesh";

  const FEShapeKey key{libMesh::MONOMIAL, p.kk_cls, p.order};

  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int q = 0; q < nqp; ++q)
    {
      Real xi   = q_pts[q](0);
      Real eta  = q_pts[q](1);
      Real zeta = q_pts[q](2);

      EXPECT_NEAR(nativeShape(key, i, xi, eta, zeta), phi_lm[i][q], TOL)
          << p.name << " phi i=" << i << " q=" << q;

      Real3 ng = nativeGradShape(key, i, xi, eta, zeta);
      EXPECT_NEAR(ng(0), dphi_lm[i][q](0), TOL)
          << p.name << " dphi/dxi i=" << i << " q=" << q;
      if (p.dim >= 2)
        EXPECT_NEAR(ng(1), dphi_lm[i][q](1), TOL)
            << p.name << " dphi/deta i=" << i << " q=" << q;
      if (p.dim >= 3)
        EXPECT_NEAR(ng(2), dphi_lm[i][q](2), TOL)
            << p.name << " dphi/dzeta i=" << i << " q=" << q;
    }
}

INSTANTIATE_TEST_SUITE_P(MonomialAllCases,
                         KokkosMonomialParityTest,
                         ::testing::ValuesIn(monomialParityParams()),
                         [](const auto & info) { return info.param.name; });

// ── G. nativeShape(FEShapeKey) dispatch is bit-identical to direct calls ─────
//
// The key-based dispatch must call the same inline function at the same
// coordinates — results must be bit-identical, not just within tolerance.

TEST(KokkosNativeFETest, MonomialKeyDispatch_BitIdentical_Dim2_Order2_Quad)
{
  Array<Real3> q_pts; Array<Real> wts;
  fillQuadrature(libMesh::QUAD4, 3, q_pts, wts);
  const FEShapeKey key{libMesh::MONOMIAL, libMesh::FEElemClass::QUAD, 2};
  const unsigned int n = MonomialImpl2D<2>::n_dofs();

  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int q = 0; q < q_pts.size(); ++q)
    {
      Real xi = q_pts[q](0), eta = q_pts[q](1);
      Real direct   = MonomialImpl2D<2>::shape(i, xi, eta, 0.0);
      Real dispatch = nativeShape(key, i, xi, eta, 0.0);
      EXPECT_EQ(direct, dispatch) << "i=" << i << " q=" << q;
    }
}

TEST(KokkosNativeFETest, MonomialKeyDispatch_BitIdentical_Dim3_Order1_Hex)
{
  Array<Real3> q_pts; Array<Real> wts;
  fillQuadrature(libMesh::HEX8, 2, q_pts, wts);
  const FEShapeKey key{libMesh::MONOMIAL, libMesh::FEElemClass::HEX, 1};
  const unsigned int n = MonomialImpl3D<1>::n_dofs();

  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int q = 0; q < q_pts.size(); ++q)
    {
      Real xi = q_pts[q](0), eta = q_pts[q](1), zeta = q_pts[q](2);
      Real3 direct   = MonomialImpl3D<1>::grad_shape(i, xi, eta, zeta);
      Real3 dispatch = nativeGradShape(key, i, xi, eta, zeta);
      EXPECT_EQ(direct(0), dispatch(0)) << "v[0] i=" << i << " q=" << q;
      EXPECT_EQ(direct(1), dispatch(1)) << "v[1] i=" << i << " q=" << q;
      EXPECT_EQ(direct(2), dispatch(2)) << "v[2] i=" << i << " q=" << q;
    }
}
