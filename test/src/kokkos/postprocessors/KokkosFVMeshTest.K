//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosFVMeshTest.h"

#include "FaceInfo.h"
#include "ElemInfo.h"
#include "FEProblemBase.h"
#include "MooseMesh.h"

#include <algorithm>

registerMooseObject("MooseTestApp", KokkosFVMeshTest);

namespace
{
bool
equal(const Moose::Kokkos::Real3 & a, const Point & b)
{
  return a(0) == b(0) && a(1) == b(1) && a(2) == b(2);
}
}

InputParameters
KokkosFVMeshTest::validParams()
{
  auto params = GeneralPostprocessor::validParams();
  params.addClassDescription("Checks the Kokkos finite volume mesh data against the MOOSE face and "
                             "element info, and returns the number of failed checks.");
  // Summing six area-weighted unit normals of order-one magnitude leaves round-off well below this
  params.addParam<Real>("closure_tolerance",
                        1e-10,
                        "Relative tolerance of the sum of the outward area-weighted normals of a "
                        "cell. Only meaningful in Cartesian coordinates.");
  return params;
}

KokkosFVMeshTest::KokkosFVMeshTest(const InputParameters & parameters)
  : GeneralPostprocessor(parameters), _closure_tolerance(getParam<Real>("closure_tolerance"))
{
}

void
KokkosFVMeshTest::execute()
{
  checkHost();
  checkDevice();
}

void
KokkosFVMeshTest::checkHost()
{
  const auto & fv_mesh = _fe_problem.getKokkosFVMesh();
  const auto & kokkos_mesh = *_fe_problem.mesh().getKokkosMesh();
  const auto & face_info = _fe_problem.mesh().faceInfo();

  auto check = [this](const bool ok, const std::string & what)
  {
    if (!ok)
    {
      ++_failures;
      _console << "KokkosFVMeshTest failure: " << what << std::endl;
    }
  };

  check(fv_mesh.getNumFaces() == face_info.size(), "number of faces");

  for (const auto * const fi : face_info)
  {
    const auto face = fv_mesh.getContiguousFaceID(*fi);
    const auto id = std::to_string(fi->id());
    check(&fv_mesh.getFaceInfo(face) == fi, "face round trip of FaceInfo " + id);

    const auto elem = kokkos_mesh.getContiguousElementID(fi->elemPtr());
    check(fv_mesh.getElem(face) == elem, "elem of FaceInfo " + id);
    check(fv_mesh.getElemSide(face) == fi->elemSideID(), "elem side of FaceInfo " + id);
    check(fv_mesh.getElemVolume(elem) == fi->elemInfo()->volume() * fi->elemInfo()->coordFactor(),
          "elem volume of FaceInfo " + id);
    check(equal(fv_mesh.getElemCentroid(elem), fi->elemCentroid()),
          "elem centroid of FaceInfo " + id);

    if (const auto * const neighbor = fi->neighborPtr())
    {
      const auto neighbor_elem = kokkos_mesh.getContiguousElementID(neighbor);
      check(fv_mesh.hasNeighbor(face), "neighbor existence of FaceInfo " + id);
      check(fv_mesh.getNeighbor(face) == neighbor_elem, "neighbor of FaceInfo " + id);
      check(fv_mesh.getNeighborSide(face) == fi->neighborSideID(),
            "neighbor side of FaceInfo " + id);
      check(fv_mesh.getElemVolume(neighbor_elem) ==
                fi->neighborInfo()->volume() * fi->neighborInfo()->coordFactor(),
            "neighbor volume of FaceInfo " + id);
    }
    else
      check(!fv_mesh.hasNeighbor(face), "missing neighbor of FaceInfo " + id);

    check(fv_mesh.getArea(face) == fi->faceArea() * fi->faceCoord(), "area of FaceInfo " + id);
    check(equal(fv_mesh.getNormal(face), fi->normal()),
          "normal of FaceInfo " + id);
    check(equal(fv_mesh.getCentroid(face), fi->faceCentroid()),
          "centroid of FaceInfo " + id);
    check(fv_mesh.getGC(face) == fi->gC(), "gC of FaceInfo " + id);
    check(fv_mesh.getDCNMag(face) == fi->dCNMag(), "dCN magnitude of FaceInfo " + id);
  }

  // The ranges cover the faces in order, and every face has the properties of its range
  const auto pid = processor_id();
  ContiguousFaceID next_face = 0;
  const auto & ranges = fv_mesh.getFaceRanges();
  for (const auto r : index_range(ranges))
  {
    const auto & range = ranges[r];
    const auto id = std::to_string(r);
    check(range.begin == next_face && range.end > range.begin, "extent of range " + id);
    next_face = range.end;

    for (const auto face : make_range(range.begin, range.end))
    {
      const auto & fi = fv_mesh.getFaceInfo(face);
      const auto * const neighbor = fi.neighborPtr();
      using Moose::Kokkos::FaceKind;
      const auto kind = !neighbor ? FaceKind::BOUNDARY
                        : neighbor->subdomain_id() != fi.elemSubdomainID() ? FaceKind::INTERFACE
                        : (fi.elem().processor_id() == pid && neighbor->processor_id() == pid)
                            ? FaceKind::INTERIOR
                            : FaceKind::PROCESSOR;
      check(range.kind == kind, "kind of a face of range " + id);
      check(range.elem_subdomain == fi.elemSubdomainID(), "elem block of a face of range " + id);
      check(range.neighbor_subdomain ==
                (neighbor ? neighbor->subdomain_id() : Moose::INVALID_BLOCK_ID),
            "neighbor block of a face of range " + id);
      check(range.elem_owned == (fi.elem().processor_id() == pid),
            "elem ownership of a face of range " + id);
      check(range.neighbor_owned == (neighbor && neighbor->processor_id() == pid),
            "neighbor ownership of a face of range " + id);
      if (kind == FaceKind::BOUNDARY || kind == FaceKind::INTERFACE)
      {
        auto sorted_ids = [](std::vector<BoundaryID> ids)
        {
          std::sort(ids.begin(), ids.end());
          return ids;
        };
        const auto & mesh = _fe_problem.mesh();
        check(fv_mesh.getFaceRangeBoundaryIDs(r, false) ==
                  sorted_ids(mesh.getBoundaryIDs(fi.elemPtr(), fi.elemSideID())),
              "elem side boundary IDs of a face of range " + id);
        if (neighbor)
          check(fv_mesh.getFaceRangeBoundaryIDs(r, true) ==
                    sorted_ids(mesh.getBoundaryIDs(neighbor, fi.neighborSideID())),
                "neighbor side boundary IDs of a face of range " + id);
      }
    }
  }
  check(next_face == fv_mesh.getNumFaces(), "ranges cover all faces");

  for (const auto elem : _fe_problem.mesh().getMesh().active_local_element_ptr_range())
  {
    const auto contiguous_elem = kokkos_mesh.getContiguousElementID(elem);
    for (const auto side : elem->side_index_range())
    {
      const auto * const fi = _fe_problem.mesh().faceInfo(elem, side);
      const auto face = fv_mesh.getElemFace(contiguous_elem, side);
      // The face of a side is stored by the cell that owns its FaceInfo, the elem cell, or else
      // by the neighbor
      const auto * const expected =
          fi ? fi : _fe_problem.mesh().faceInfo(elem->neighbor_ptr(side),
                                                elem->neighbor_ptr(side)->which_neighbor_am_i(elem));
      check(&fv_mesh.getFaceInfo(face) == expected,
            "face of side " + std::to_string(side) + " of element " + std::to_string(elem->id()));
    }
    for (const auto side : make_range(elem->n_sides(), fv_mesh.getMaxSidesPerElem()))
      check(fv_mesh.getElemFace(contiguous_elem, side) == libMesh::DofObject::invalid_id,
            "padding side " + std::to_string(side) + " of element " + std::to_string(elem->id()));
  }
}

void
KokkosFVMeshTest::checkDevice()
{
  const auto fv_mesh = _fe_problem.getKokkosFVMesh();
  const auto max_sides = fv_mesh.getMaxSidesPerElem();
  const auto tolerance = _closure_tolerance;

  ContiguousElementID open_cells = 0;
  ::Kokkos::parallel_reduce(
      ::Kokkos::RangePolicy<ExecSpace>(0, fv_mesh.getNumLocalElements()),
      KOKKOS_LAMBDA(const ContiguousElementID elem, ContiguousElementID & open) {
        Moose::Kokkos::Real3 closure;
        Real area = 0;
        for (unsigned int side = 0; side < max_sides; ++side)
        {
          const auto face = fv_mesh.getElemFace(elem, side);
          if (face == libMesh::DofObject::invalid_id)
            continue;
          const Real sign = fv_mesh.getElem(face) == elem ? 1 : -1;
          closure += sign * fv_mesh.getArea(face) * fv_mesh.getNormal(face);
          area += fv_mesh.getArea(face);
        }
        if (closure.norm() > tolerance * area)
          ++open;
      },
      open_cells);

  if (open_cells)
  {
    _failures += open_cells;
    _console << "KokkosFVMeshTest failure: " << open_cells << " cells whose faces do not close"
             << std::endl;
  }
}
