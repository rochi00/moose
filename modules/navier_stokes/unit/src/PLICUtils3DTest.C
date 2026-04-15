//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "PLICUtils3D.h"
#include "AppFactory.h"
#include "Factory.h"
#include "GeneratedMeshGenerator.h"
#include "MeshGeneratorMesh.h"
#include "MooseMain.h"

#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"

namespace
{

std::vector<libMesh::Point>
unitCubeVertices()
{
  using libMesh::Point;
  return {
      Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0),
      Point(0, 0, 1), Point(1, 0, 1), Point(1, 1, 1), Point(0, 1, 1)};
}

std::vector<std::vector<unsigned int>>
unitCubeFaces()
{
  return {
      {0, 3, 2, 1}, {0, 1, 5, 4}, {1, 2, 6, 5},
      {2, 3, 7, 6}, {0, 4, 7, 3}, {4, 5, 6, 7}};
}

/**
 * Helper struct that keeps app + mesh alive for the duration of a test.
 * Must be destroyed (go out of scope) before libMesh shuts down —
 * so it must NOT be static.
 */
struct HexElemFixture
{
  std::shared_ptr<MooseApp> app;
  std::shared_ptr<MooseMesh> mesh;
  const libMesh::Elem * elem;

  HexElemFixture()
  {
    const char * argv[2] = {"foo", "\0"};
    app = Moose::createMooseApp("NavierStokesUnitApp", 1, (char **)argv);
    auto * factory = &app->getFactory();

    InputParameters mesh_params = factory->getValidParams("MeshGeneratorMesh");
    mesh = factory->create<MeshGeneratorMesh>("MeshGeneratorMesh", "moose_mesh", mesh_params);
    app->actionWarehouse().mesh() = mesh;

    InputParameters gen_params = factory->getValidParams("GeneratedMeshGenerator");
    gen_params.set<unsigned int>("nx") = 1;
    gen_params.set<unsigned int>("ny") = 1;
    gen_params.set<unsigned int>("nz") = 1;
    gen_params.set<MooseEnum>("dim") = "3";
    auto mesh_gen = factory->create<GeneratedMeshGenerator>(
        "GeneratedMeshGenerator", "mesh_gen", gen_params);

    auto lm_mesh = mesh_gen->generate();
    mesh->setMeshBase(std::move(lm_mesh));
    mesh->prepare(nullptr);

    elem = *mesh->getMesh().active_local_element_ptr_range().begin();
  }
};

} // namespace

TEST(PLICUtils3DTest, clippedFaceAreaHalfSquare)
{
  using libMesh::Point;
  using libMesh::VectorValue;

  const std::vector<Point> face = {
      Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)};
  const VectorValue<Real> normal(1, 0, 0);

  EXPECT_NEAR(NS::PLIC::clippedFaceArea(face, normal, 0.5), 0.5, 1e-14);
}

TEST(PLICUtils3DTest, clippedVolumeHalfCube)
{
  using libMesh::VectorValue;

  const VectorValue<Real> normal(1, 0, 0);
  EXPECT_NEAR(
      NS::PLIC::clippedVolume(unitCubeVertices(), unitCubeFaces(), normal, 0.5), 0.5, 1e-14);
}

TEST(PLICUtils3DTest, computePlaneOffset3DRecoversTargetAlpha)
{
  using libMesh::VectorValue;

  VectorValue<Real> normal(1, 1, 1);
  normal /= normal.norm();

  const Real alpha = 0.37;
  const Real offset =
      NS::PLIC::computePlaneOffset3D(unitCubeVertices(), unitCubeFaces(), 1.0, normal, alpha);
  const Real reconstructed = NS::PLIC::truncatedVolumeFraction3D(
      unitCubeVertices(), unitCubeFaces(), 1.0, normal, offset);

  EXPECT_NEAR(reconstructed, alpha, 1e-12);
}

TEST(PLICUtils3DTest, elemGeometryExtractsHexFaces)
{
  using libMesh::VectorValue;

  // Fixture owns app+mesh and is destroyed at end of this scope,
  // before libMesh shuts down — no leak.
  HexElemFixture fixture;

  std::vector<libMesh::Point> vertices;
  std::vector<std::vector<unsigned int>> faces;
  NS::PLIC::elemGeometry(*fixture.elem, vertices, faces);

  ASSERT_EQ(vertices.size(), 8);
  ASSERT_EQ(faces.size(), 6);

  const VectorValue<Real> normal(1, 0, 0);
  EXPECT_NEAR(NS::PLIC::truncatedVolumeFraction3D(vertices, faces, 1.0, normal, 0.5), 0.5, 1e-14);
}

TEST(PLICUtils3DTest, timeIntegratedSubmergedAreaMatchesExactIntegral)
{
  using libMesh::Point;
  using libMesh::VectorValue;

  const std::vector<Point> face = {
      Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0)};
  const VectorValue<Real> normal(1, 0, 0);

  EXPECT_NEAR(NS::PLIC::timeIntegratedSubmergedArea(face, normal, 0.0, 0.5, 1.0, 1.0),
              0.5,
              1e-14);
}
