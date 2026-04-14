//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest_include.h"

#include "libmesh/gpu/kokkos_fe_types.h"

#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"

using namespace libMesh::Kokkos;

// ── getSideTopology ───────────────────────────────────────────────────────────

TEST(KokkosFETypesTest, GetSideTopology_1D)
{
  // 1D edge "sides" are vertex nodes; represented as EDGE2
  EXPECT_EQ(getSideTopology(libMesh::EDGE2), libMesh::EDGE2);
  EXPECT_EQ(getSideTopology(libMesh::EDGE3), libMesh::EDGE2);
}

TEST(KokkosFETypesTest, GetSideTopology_2D_FirstOrder)
{
  // First-order 2D elements have linear edge sides
  EXPECT_EQ(getSideTopology(libMesh::TRI3),  libMesh::EDGE2);
  EXPECT_EQ(getSideTopology(libMesh::QUAD4), libMesh::EDGE2);
}

TEST(KokkosFETypesTest, GetSideTopology_2D_SecondOrder)
{
  // Second-order 2D elements have quadratic edge sides
  EXPECT_EQ(getSideTopology(libMesh::TRI6),  libMesh::EDGE3);
  EXPECT_EQ(getSideTopology(libMesh::QUAD8), libMesh::EDGE3);
  EXPECT_EQ(getSideTopology(libMesh::QUAD9), libMesh::EDGE3);
}

TEST(KokkosFETypesTest, GetSideTopology_3D_FirstOrder)
{
  EXPECT_EQ(getSideTopology(libMesh::TET4), libMesh::TRI3);
  EXPECT_EQ(getSideTopology(libMesh::HEX8), libMesh::QUAD4);
}

TEST(KokkosFETypesTest, GetSideTopology_3D_SecondOrder)
{
  EXPECT_EQ(getSideTopology(libMesh::TET10), libMesh::TRI6);
  EXPECT_EQ(getSideTopology(libMesh::HEX20), libMesh::QUAD8);
  EXPECT_EQ(getSideTopology(libMesh::HEX27), libMesh::QUAD9);
}

// ── classFromTopology ─────────────────────────────────────────────────────────

TEST(KokkosFETypesTest, ClassFromTopology_1D)
{
  EXPECT_EQ(classFromTopology(libMesh::EDGE2), libMesh::FEElemClass::EDGE);
  EXPECT_EQ(classFromTopology(libMesh::EDGE3), libMesh::FEElemClass::EDGE);
}

TEST(KokkosFETypesTest, ClassFromTopology_2D)
{
  EXPECT_EQ(classFromTopology(libMesh::TRI3),  libMesh::FEElemClass::TRI);
  EXPECT_EQ(classFromTopology(libMesh::TRI6),  libMesh::FEElemClass::TRI);
  EXPECT_EQ(classFromTopology(libMesh::QUAD4), libMesh::FEElemClass::QUAD);
  EXPECT_EQ(classFromTopology(libMesh::QUAD8), libMesh::FEElemClass::QUAD);
  EXPECT_EQ(classFromTopology(libMesh::QUAD9), libMesh::FEElemClass::QUAD);
}

TEST(KokkosFETypesTest, ClassFromTopology_3D)
{
  EXPECT_EQ(classFromTopology(libMesh::TET4),  libMesh::FEElemClass::TET);
  EXPECT_EQ(classFromTopology(libMesh::TET10), libMesh::FEElemClass::TET);
  EXPECT_EQ(classFromTopology(libMesh::HEX8),  libMesh::FEElemClass::HEX);
  EXPECT_EQ(classFromTopology(libMesh::HEX20), libMesh::FEElemClass::HEX);
  EXPECT_EQ(classFromTopology(libMesh::HEX27), libMesh::FEElemClass::HEX);
}

// ── nDofs(FEShapeKey) — LAGRANGE ─────────────────────────────────────────────

TEST(KokkosFETypesTest, NDofs_FEShapeKey_Lagrange_1D)
{
  EXPECT_EQ(nDofs(FEShapeKey{libMesh::LAGRANGE, libMesh::FEElemClass::EDGE, 1}), 2u);
  EXPECT_EQ(nDofs(FEShapeKey{libMesh::LAGRANGE, libMesh::FEElemClass::EDGE, 2}), 3u);
}

TEST(KokkosFETypesTest, NDofs_FEShapeKey_Lagrange_2D)
{
  EXPECT_EQ(nDofs(FEShapeKey{libMesh::LAGRANGE, libMesh::FEElemClass::TRI,  1}),  3u);
  EXPECT_EQ(nDofs(FEShapeKey{libMesh::LAGRANGE, libMesh::FEElemClass::TRI,  2}),  6u);
  EXPECT_EQ(nDofs(FEShapeKey{libMesh::LAGRANGE, libMesh::FEElemClass::QUAD, 1}),  4u);
  EXPECT_EQ(nDofs(FEShapeKey{libMesh::LAGRANGE, libMesh::FEElemClass::QUAD, 2}),  9u);
}

TEST(KokkosFETypesTest, NDofs_FEShapeKey_Lagrange_3D)
{
  EXPECT_EQ(nDofs(FEShapeKey{libMesh::LAGRANGE, libMesh::FEElemClass::TET, 1}),  4u);
  EXPECT_EQ(nDofs(FEShapeKey{libMesh::LAGRANGE, libMesh::FEElemClass::TET, 2}), 10u);
  EXPECT_EQ(nDofs(FEShapeKey{libMesh::LAGRANGE, libMesh::FEElemClass::HEX, 1}),  8u);
  EXPECT_EQ(nDofs(FEShapeKey{libMesh::LAGRANGE, libMesh::FEElemClass::HEX, 2}), 27u);
}

// ── nDofs(FEShapeKey) — MONOMIAL ─────────────────────────────────────────────
// 1D: p+1  |  2D: (p+1)(p+2)/2  |  3D: (p+1)(p+2)(p+3)/6

TEST(KokkosFETypesTest, NDofs_FEShapeKey_Monomial_1D)
{
  for (unsigned int p = 0; p <= 5; ++p)
    EXPECT_EQ(nDofs(FEShapeKey{libMesh::MONOMIAL, libMesh::FEElemClass::EDGE, p}), p + 1) << "p=" << p;
}

TEST(KokkosFETypesTest, NDofs_FEShapeKey_Monomial_2D_TRI)
{
  for (unsigned int p = 0; p <= 5; ++p)
    EXPECT_EQ(nDofs(FEShapeKey{libMesh::MONOMIAL, libMesh::FEElemClass::TRI, p}),
              (p + 1) * (p + 2) / 2) << "p=" << p;
}

TEST(KokkosFETypesTest, NDofs_FEShapeKey_Monomial_2D_QUAD)
{
  // TRI and QUAD share the same Dim2 formula (total-degree polynomial space)
  for (unsigned int p = 0; p <= 5; ++p)
    EXPECT_EQ(nDofs(FEShapeKey{libMesh::MONOMIAL, libMesh::FEElemClass::QUAD, p}),
              (p + 1) * (p + 2) / 2) << "p=" << p;
}

TEST(KokkosFETypesTest, NDofs_FEShapeKey_Monomial_3D_TET)
{
  for (unsigned int p = 0; p <= 5; ++p)
    EXPECT_EQ(nDofs(FEShapeKey{libMesh::MONOMIAL, libMesh::FEElemClass::TET, p}),
              (p + 1) * (p + 2) * (p + 3) / 6) << "p=" << p;
}

TEST(KokkosFETypesTest, NDofs_FEShapeKey_Monomial_3D_HEX)
{
  // TET and HEX share the same Dim3 formula
  for (unsigned int p = 0; p <= 5; ++p)
    EXPECT_EQ(nDofs(FEShapeKey{libMesh::MONOMIAL, libMesh::FEElemClass::HEX, p}),
              (p + 1) * (p + 2) * (p + 3) / 6) << "p=" << p;
}
