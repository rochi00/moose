//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosFVMesh.h"

#include "FaceInfo.h"
#include "ElemInfo.h"
#include "MooseMesh.h"

#include <algorithm>
#include <tuple>

namespace Moose::Kokkos
{

void
FVMesh::update()
{
  const auto & kokkos_mesh = *_mesh.getKokkosMesh();
  mooseAssert(kokkos_mesh.initialized(),
              "The Kokkos mesh must be initialized before the finite volume mesh data");

  const auto & face_info = _mesh.faceInfo();

  _num_faces = face_info.size();
  _num_local_elems = kokkos_mesh.getNumLocalElements();
  _max_sides = _mesh.getMaxSidesPerElem();
  const auto num_elems = kokkos_mesh.getNumLocalAndPossiblyOneNeighborLayerGhostElements();
  _num_elems = num_elems;

  _host_maps = std::make_shared<HostMaps>();
  auto & maps = *_host_maps;

  // Order the faces into ranges whose faces share every property a face loop branches on, and
  // within a range by their owned cell, so that consecutive threads touch nearby cells
  const auto pid = _mesh.processor_id();
  struct FaceKey
  {
    std::tuple<SubdomainID,
               FaceKind,
               SubdomainID,
               bool,
               bool,
               std::vector<BoundaryID>,
               std::vector<BoundaryID>>
        range;
    ContiguousElementID owned_elem;
    unsigned int side;
    const FaceInfo * fi;
  };
  std::vector<FaceKey> keys;
  keys.reserve(face_info.size());
  for (const auto * const fi : face_info)
  {
    const auto * const neighbor = fi->neighborPtr();
    const bool elem_owned = fi->elem().processor_id() == pid;
    const bool neighbor_owned = neighbor && neighbor->processor_id() == pid;
    const auto neighbor_subdomain = neighbor ? neighbor->subdomain_id() : Moose::INVALID_BLOCK_ID;

    FaceKind kind;
    if (!neighbor)
      kind = FaceKind::BOUNDARY;
    else if (neighbor->subdomain_id() != fi->elemSubdomainID())
      kind = FaceKind::INTERFACE;
    else if (elem_owned && neighbor_owned)
      kind = FaceKind::INTERIOR;
    else
      kind = FaceKind::PROCESSOR;

    // A variable can live on one side of a face only at a block interface or a mesh boundary,
    // which is where boundary conditions apply. A sideset holds element sides, so the IDs are
    // kept per side.
    std::vector<BoundaryID> elem_boundary_ids, neighbor_boundary_ids;
    if (kind == FaceKind::BOUNDARY || kind == FaceKind::INTERFACE)
    {
      elem_boundary_ids = _mesh.getBoundaryIDs(fi->elemPtr(), fi->elemSideID());
      std::sort(elem_boundary_ids.begin(), elem_boundary_ids.end());
      if (neighbor)
      {
        neighbor_boundary_ids = _mesh.getBoundaryIDs(neighbor, fi->neighborSideID());
        std::sort(neighbor_boundary_ids.begin(), neighbor_boundary_ids.end());
      }
    }

    const auto * const owned = elem_owned ? fi->elemPtr() : neighbor;
    keys.push_back({{fi->elemSubdomainID(),
                     kind,
                     neighbor_subdomain,
                     elem_owned,
                     neighbor_owned,
                     std::move(elem_boundary_ids),
                     std::move(neighbor_boundary_ids)},
                    kokkos_mesh.getContiguousElementID(owned),
                    elem_owned ? fi->elemSideID() : fi->neighborSideID(),
                    fi});
  }
  std::sort(keys.begin(),
            keys.end(),
            [](const FaceKey & a, const FaceKey & b)
            {
              return std::tie(a.range, a.owned_elem, a.side) <
                     std::tie(b.range, b.owned_elem, b.side);
            });

  maps.face_infos.resize(keys.size());
  for (const auto face : index_range(keys))
  {
    maps.face_infos[face] = keys[face].fi;
    if (!face || keys[face].range != keys[face - 1].range)
    {
      const auto & [elem_subdomain,
                    kind,
                    neighbor_subdomain,
                    elem_owned,
                    neighbor_owned,
                    elem_ids,
                    neighbor_ids] = keys[face].range;
      maps.face_ranges.push_back({face,
                                  face,
                                  kind,
                                  elem_subdomain,
                                  neighbor_subdomain,
                                  elem_owned,
                                  neighbor_owned});
      maps.face_range_elem_boundary_ids.push_back(elem_ids);
      maps.face_range_neighbor_boundary_ids.push_back(neighbor_ids);
    }
    maps.face_ranges.back().end = face + 1;
  }

  dof_id_type max_face_info_id = 0;
  for (const auto * const fi : face_info)
    max_face_info_id = std::max(max_face_info_id, fi->id());
  maps.face_info_to_face.assign(max_face_info_id + 1, libMesh::DofObject::invalid_id);

  _face_elem.create(_num_faces);
  _face_neighbor.create(_num_faces);
  _face_elem_side.create(_num_faces);
  _face_neighbor_side.create(_num_faces);
  _face_area.create(_num_faces);
  _face_normal.create(_num_faces);
  _face_centroid.create(_num_faces);
  _face_gc.create(_num_faces);
  _face_dcn_mag.create(_num_faces);

  _elem_volume.create(num_elems);
  _elem_centroid.create(num_elems);
  _elem_volume = 0;
  _elem_centroid = Real3();

  _elem_faces.create(_num_local_elems, _max_sides);
  _elem_faces = libMesh::DofObject::invalid_id;

  // Record the volume and centroid of a cell, and the face on its side if the cell is owned
  auto add_cell = [&](const ElemInfo & elem_info,
                      const ContiguousElementID elem,
                      const unsigned int side,
                      const ContiguousFaceID face)
  {
    _elem_volume[elem] = elem_info.volume() * elem_info.coordFactor();
    _elem_centroid[elem] = elem_info.centroid();

    if (elem >= _num_local_elems)
      return;

    auto & entry = _elem_faces(elem, side);
    // A cell side facing several finer cells has several faces, which the cell-to-face table
    // cannot hold
    if (entry != libMesh::DofObject::invalid_id)
      mooseError("The Kokkos linear finite volume objects do not support meshes where a cell side "
                 "has more than one face, such as adaptively refined meshes.");
    entry = face;
  };

  for (const auto face : make_range(_num_faces))
  {
    const auto & fi = *maps.face_infos[face];
    maps.face_info_to_face[fi.id()] = face;

    const auto elem = kokkos_mesh.getContiguousElementID(fi.elemPtr());
    _face_elem[face] = elem;
    _face_elem_side[face] = fi.elemSideID();
    add_cell(*fi.elemInfo(), elem, fi.elemSideID(), face);

    if (const auto * const neighbor = fi.neighborPtr())
    {
      const auto neighbor_elem = kokkos_mesh.getContiguousElementID(neighbor);
      _face_neighbor[face] = neighbor_elem;
      _face_neighbor_side[face] = fi.neighborSideID();
      add_cell(*fi.neighborInfo(), neighbor_elem, fi.neighborSideID(), face);
    }
    else
    {
      _face_neighbor[face] = libMesh::DofObject::invalid_id;
      _face_neighbor_side[face] = libMesh::invalid_uint;
    }

    _face_area[face] = fi.faceArea() * fi.faceCoord();
    _face_normal[face] = fi.normal();
    _face_centroid[face] = fi.faceCentroid();
    _face_gc[face] = fi.gC();
    _face_dcn_mag[face] = fi.dCNMag();
  }

#ifndef NDEBUG
  for (const auto elem : _mesh.getMesh().active_local_element_ptr_range())
  {
    const auto contiguous_elem = kokkos_mesh.getContiguousElementID(elem);
    for (const auto side : elem->side_index_range())
      mooseAssert(_elem_faces(contiguous_elem, side) != libMesh::DofObject::invalid_id,
                  "Every side of an owned cell must have a face");
  }
#endif

  _face_elem.copyToDevice();
  _face_neighbor.copyToDevice();
  _face_elem_side.copyToDevice();
  _face_neighbor_side.copyToDevice();
  _face_area.copyToDevice();
  _face_normal.copyToDevice();
  _face_centroid.copyToDevice();
  _face_gc.copyToDevice();
  _face_dcn_mag.copyToDevice();
  _elem_volume.copyToDevice();
  _elem_centroid.copyToDevice();
  _elem_faces.copyToDevice();

  _initialized = true;

  if (_needs_side_data)
    buildSideData();
}

void
FVMesh::requestSideData() const
{
  mooseAssert(_initialized, "The finite volume mesh data must be built first");

  if (_needs_side_data)
    return;
  _needs_side_data = true;
  buildSideData();
}

void
FVMesh::buildSideData() const
{
  _side_neighbor.create(_num_local_elems, _max_sides);
  _side_weight.create(_num_local_elems, _max_sides);
  _side_area_normal.create(_num_local_elems, _max_sides);
  _side_neighbor = libMesh::DofObject::invalid_id;
  _side_weight = 0;
  _side_area_normal = Real3();

  for (const auto face : make_range(_num_faces))
  {
    const auto elem = _face_elem[face];
    const auto neighbor = _face_neighbor[face];
    const auto area_normal = _face_area[face] * _face_normal[face];

    if (elem < _num_local_elems)
    {
      const auto side = _face_elem_side[face];
      _side_neighbor(elem, side) = neighbor;
      _side_weight(elem, side) = _face_gc[face];
      _side_area_normal(elem, side) = area_normal;
    }
    if (neighbor != libMesh::DofObject::invalid_id && neighbor < _num_local_elems)
    {
      const auto side = _face_neighbor_side[face];
      _side_neighbor(neighbor, side) = elem;
      _side_weight(neighbor, side) = 1 - _face_gc[face];
      _side_area_normal(neighbor, side) = -area_normal;
    }
  }

  _side_neighbor.copyToDevice();
  _side_weight.copyToDevice();
  _side_area_normal.copyToDevice();
}

ContiguousFaceID
FVMesh::getContiguousFaceID(const FaceInfo & fi) const
{
  const auto & face_info_to_face = _host_maps->face_info_to_face;
  mooseAssert(fi.id() < face_info_to_face.size() &&
                  face_info_to_face[fi.id()] != libMesh::DofObject::invalid_id,
              "The face has no owned cell on either side");

  return face_info_to_face[fi.id()];
}

void
FaceRangeList::create(const std::vector<FaceRange> & ranges)
{
  std::vector<ContiguousFaceID> begin;
  std::vector<MOOSE_KOKKOS_INDEX_TYPE> offset{0};
  for (const auto & range : ranges)
  {
    begin.push_back(range.begin);
    offset.push_back(offset.back() +
                     libMesh::cast_int<MOOSE_KOKKOS_INDEX_TYPE>(range.end - range.begin));
  }

  _size = offset.back();
  _num_ranges = ranges.size();
  _begin = begin;
  _offset = offset;
}

MOOSE_KOKKOS_INDEX_TYPE
FaceRangeList::offset(const unsigned int range) const
{
  return _offset[range];
}

ContiguousFaceID
FaceRangeList::begin(const unsigned int range) const
{
  return _begin[range];
}

} // namespace Moose::Kokkos
