//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosMesh.h"

#include <memory>
#include <set>
#include <vector>

class FaceInfo;

using ContiguousFaceID = dof_id_type;

namespace Moose::Kokkos
{

/// Kind of a face, from the cells on its two sides
enum class FaceKind : unsigned int
{
  /// Both cells owned, in the same block
  INTERIOR,
  /// One cell owned and one ghost, in the same block
  PROCESSOR,
  /// Cells in different blocks
  INTERFACE,
  /// No neighbor cell
  BOUNDARY
};

/**
 * A contiguous range of faces that share every property a face loop branches on: the blocks of
 * their cells, their kind, which of their cells are owned, and the boundary IDs of each of their
 * sides
 */
struct FaceRange
{
  /// First face of the range
  ContiguousFaceID begin;
  /// One past the last face of the range
  ContiguousFaceID end;
  /// Kind of the faces
  FaceKind kind;
  /// Block of the elem cells
  SubdomainID elem_subdomain;
  /// Block of the neighbor cells, Moose::INVALID_BLOCK_ID on a mesh boundary
  SubdomainID neighbor_subdomain;
  /// Whether the elem cells are owned
  bool elem_owned;
  /// Whether the neighbor cells are owned
  bool neighbor_owned;
};

/**
 * Face and cell data of the Kokkos linear finite volume objects, built from the MOOSE FaceInfo and
 * ElemInfo objects so that the geometry is identical to the one of the host linear finite volume
 * objects.
 *
 * The faces are every face of MooseMesh::faceInfo(), those with an owned cell on at least one side.
 * They are numbered in ranges (see FaceRange) and, within a range, in the order of their owned
 * cells. Each face stores the contiguous IDs (see Mesh) of its two cells, so a
 * cell can be a one-layer ghost. The orientation is the one of FaceInfo: the normal points out of
 * the elem cell, and the neighbor cell is invalid on mesh boundaries. Cells store the faces around
 * them in a cell-fastest table, so that consecutive threads of a cell loop read consecutive
 * entries.
 *
 * The data is kept on both host and device, so that setup code can read it on host.
 */
class FVMesh
{
public:
  /**
   * Constructor
   * @param mesh The MOOSE mesh
   */
  FVMesh(const MooseMesh & mesh) : _mesh(mesh) {}

  /**
   * Rebuild the face and cell data. The FaceInfo and ElemInfo objects of the MOOSE mesh and the
   * Kokkos mesh must be up to date.
   */
  void update();

  /**
   * Get whether the data was built
   */
  [[nodiscard]] bool initialized() const { return _initialized; }

  /**
   * Get the number of faces
   */
  ContiguousFaceID getNumFaces() const { return _num_faces; }
  /**
   * Get the number of owned cells
   */
  ContiguousElementID getNumLocalElements() const { return _num_local_elems; }
  /**
   * Get the maximum number of sides of a cell, the second dimension of the cell-to-face table
   */
  unsigned int getMaxSidesPerElem() const { return _max_sides; }

  /**
   * Build the per-side data of the owned cells (getSideNeighbor() and the following), if not built
   * yet. It is kept up to date by update() once built. Built on request, as only cell gradient
   * computations need it.
   */
  void requestSideData() const;

  /**
   * Get the contiguous ID of a face
   * @param fi The face, which must have an owned cell on at least one side
   * @returns The contiguous face ID
   */
  ContiguousFaceID getContiguousFaceID(const FaceInfo & fi) const;
  /**
   * Get the MOOSE face of a contiguous face ID (host only)
   * @param face The contiguous face ID
   * @returns The MOOSE face
   */
  const FaceInfo & getFaceInfo(const ContiguousFaceID face) const
  {
    return *_host_maps->face_infos[face];
  }
  /**
   * Get the face ranges, which cover the faces in order (host only)
   */
  const std::vector<FaceRange> & getFaceRanges() const { return _host_maps->face_ranges; }
  /**
   * Get the boundary IDs of the faces of a range on one side (host only). A sideset holds element
   * sides, so the two sides of an internal face can carry different boundary IDs. They are kept
   * only at block interfaces and mesh boundaries, where a variable can live on one side only;
   * elsewhere the sets are empty.
   * @param range The index of the range in getFaceRanges()
   * @param neighbor_side Whether to get the IDs of the neighbor cell's side instead of the elem
   * cell's
   */
  const std::vector<BoundaryID> & getFaceRangeBoundaryIDs(const unsigned int range,
                                                          const bool neighbor_side) const
  {
    return neighbor_side ? _host_maps->face_range_neighbor_boundary_ids[range]
                         : _host_maps->face_range_elem_boundary_ids[range];
  }

#ifdef MOOSE_KOKKOS_SCOPE
  /**
   * Get the contiguous ID of the elem cell of a face, out of which the face normal points
   */
  KOKKOS_FUNCTION ContiguousElementID getElem(const ContiguousFaceID face) const
  {
    return _face_elem[face];
  }
  /**
   * Get the contiguous ID of the neighbor cell of a face
   * @returns The contiguous ID, or libMesh::DofObject::invalid_id on a mesh boundary
   */
  KOKKOS_FUNCTION ContiguousElementID getNeighbor(const ContiguousFaceID face) const
  {
    return _face_neighbor[face];
  }
  /**
   * Get whether a face has a neighbor cell
   */
  KOKKOS_FUNCTION bool hasNeighbor(const ContiguousFaceID face) const
  {
    return _face_neighbor[face] != libMesh::DofObject::invalid_id;
  }
  /**
   * Get the side index of a face on its elem cell
   */
  KOKKOS_FUNCTION unsigned int getElemSide(const ContiguousFaceID face) const
  {
    return _face_elem_side[face];
  }
  /**
   * Get the side index of a face on its neighbor cell
   * @returns The side index, or libMesh::invalid_uint on a mesh boundary
   */
  KOKKOS_FUNCTION unsigned int getNeighborSide(const ContiguousFaceID face) const
  {
    return _face_neighbor_side[face];
  }
  /**
   * Get the area of a face multiplied by the coordinate transformation factor
   */
  KOKKOS_FUNCTION Real getArea(const ContiguousFaceID face) const { return _face_area[face]; }
  /**
   * Get the unit normal of a face, pointing out of its elem cell
   */
  KOKKOS_FUNCTION Real3 getNormal(const ContiguousFaceID face) const { return _face_normal[face]; }
  /**
   * Get the centroid of a face
   */
  KOKKOS_FUNCTION Real3 getCentroid(const ContiguousFaceID face) const
  {
    return _face_centroid[face];
  }
  /**
   * Get the linear interpolation weight of the elem cell of a face (FaceInfo::gC()), 1 on a mesh
   * boundary
   */
  KOKKOS_FUNCTION Real getGC(const ContiguousFaceID face) const { return _face_gc[face]; }
  /**
   * Get the distance between the cell centroids of a face, or between the elem cell centroid and
   * the face centroid on a mesh boundary (FaceInfo::dCNMag())
   */
  KOKKOS_FUNCTION Real getDCNMag(const ContiguousFaceID face) const { return _face_dcn_mag[face]; }

  /**
   * Get the volume of a cell multiplied by the coordinate transformation factor
   * @param elem The contiguous ID of an owned or ghost cell
   */
  KOKKOS_FUNCTION Real getElemVolume(const ContiguousElementID elem) const
  {
    return _elem_volume[elem];
  }
  /**
   * Get the centroid of a cell
   * @param elem The contiguous ID of an owned or ghost cell
   */
  KOKKOS_FUNCTION Real3 getElemCentroid(const ContiguousElementID elem) const
  {
    return _elem_centroid[elem];
  }
  /**
   * Get the face on a side of an owned cell
   * @param elem The contiguous ID of an owned cell
   * @param side The side index
   * @returns The contiguous face ID, or libMesh::DofObject::invalid_id past the cell's sides
   */
  KOKKOS_FUNCTION ContiguousFaceID getElemFace(const ContiguousElementID elem,
                                               const unsigned int side) const
  {
    return _elem_faces(elem, side);
  }

  /**
   * Per-side data of the owned cells, stored cell-fastest so that consecutive threads of a cell loop
   * read consecutive entries (requestSideData() must have been called)
   * @param elem The contiguous ID of an owned cell
   * @param side The side index
   */
  ///@{
  /// The cell across the side, or libMesh::DofObject::invalid_id on a mesh boundary or past the
  /// cell's sides
  KOKKOS_FUNCTION ContiguousElementID getSideNeighbor(const ContiguousElementID elem,
                                                      const unsigned int side) const
  {
    return _side_neighbor(elem, side);
  }
  /// The weight of the cell's own value in the linear interpolation to the face (FaceInfo::gC()
  /// seen from the cell); 1 on a mesh boundary
  KOKKOS_FUNCTION Real getSideWeight(const ContiguousElementID elem, const unsigned int side) const
  {
    return _side_weight(elem, side);
  }
  /// The face normal pointing out of the cell, times the coordinate-weighted face area
  KOKKOS_FUNCTION Real3 getSideAreaNormal(const ContiguousElementID elem,
                                          const unsigned int side) const
  {
    return _side_area_normal(elem, side);
  }
  ///@}
#endif

private:
  /// The MOOSE mesh
  const MooseMesh & _mesh;

  /// Whether the data was built
  bool _initialized = false;
  /// Number of faces
  ContiguousFaceID _num_faces = 0;
  /// Number of owned cells
  ContiguousElementID _num_local_elems = 0;
  /// Maximum number of sides of a cell
  unsigned int _max_sides = 0;

  /// Host-only maps between MOOSE faces and contiguous face IDs, shared by the copies of this
  /// object so that dispatch copies do not copy them
  struct HostMaps
  {
    /// MOOSE face of each contiguous face ID
    std::vector<const FaceInfo *> face_infos;
    /// Contiguous face ID by FaceInfo::id(), invalid for faces without an owned cell
    std::vector<ContiguousFaceID> face_info_to_face;
    /// Face ranges, in face order
    std::vector<FaceRange> face_ranges;
    /// Boundary IDs of the elem cell side of the faces of each range
    std::vector<std::vector<BoundaryID>> face_range_elem_boundary_ids;
    /// Boundary IDs of the neighbor cell side of the faces of each range
    std::vector<std::vector<BoundaryID>> face_range_neighbor_boundary_ids;
  };
  std::shared_ptr<HostMaps> _host_maps;

  /// Face data, indexed by contiguous face ID
  ///@{
  Array<ContiguousElementID> _face_elem;
  Array<ContiguousElementID> _face_neighbor;
  Array<unsigned int> _face_elem_side;
  Array<unsigned int> _face_neighbor_side;
  Array<Real> _face_area;
  Array<Real3> _face_normal;
  Array<Real3> _face_centroid;
  Array<Real> _face_gc;
  Array<Real> _face_dcn_mag;
  ///@}

  /// Cell data, indexed by contiguous element ID of owned and ghost cells
  ///@{
  Array<Real> _elem_volume;
  Array<Real3> _elem_centroid;
  ///@}

  /// Faces around each owned cell, indexed by (cell, side) with the cell index fastest
  Array2D<ContiguousFaceID> _elem_faces;

  /// Build the per-side data
  void buildSideData() const;

  /// Whether the per-side data was requested
  mutable bool _needs_side_data = false;
  /// Per-side data of the owned cells, indexed by (cell, side) with the cell index fastest
  ///@{
  mutable Array2D<ContiguousElementID> _side_neighbor;
  mutable Array2D<Real> _side_weight;
  mutable Array2D<Real3> _side_area_normal;
  ///@}
};

/**
 * A list of face ranges that a face loop runs over, one thread per face. Thread t handles face
 * begin(r) + t - offset(r) of the range r with offset(r) <= t < offset(r + 1).
 */
class FaceRangeList
{
public:
  /**
   * Build the list
   * @param ranges The face ranges, in the order of the threads
   */
  void create(const std::vector<FaceRange> & ranges);

  /**
   * Get the number of faces, which is the number of threads
   */
  MOOSE_KOKKOS_INDEX_TYPE size() const { return _size; }
  /**
   * Get the number of ranges
   */
  unsigned int numRanges() const { return _num_ranges; }
  /**
   * Get the first thread of a range (host)
   */
  MOOSE_KOKKOS_INDEX_TYPE offset(const unsigned int range) const;

#ifdef MOOSE_KOKKOS_SCOPE
  /**
   * Get the range of a thread
   */
  KOKKOS_FUNCTION unsigned int range(const MOOSE_KOKKOS_INDEX_TYPE tid) const
  {
    // Largest range whose first thread is not past tid; the number of ranges is small
    unsigned int low = 0, high = _num_ranges;
    while (high - low > 1)
    {
      const auto mid = (low + high) / 2;
      if (_offset[mid] <= tid)
        low = mid;
      else
        high = mid;
    }
    return low;
  }
  /**
   * Get the face of a thread
   * @param tid The thread
   * @param range The range of the thread, from range()
   */
  KOKKOS_FUNCTION ContiguousFaceID face(const MOOSE_KOKKOS_INDEX_TYPE tid,
                                        const unsigned int range) const
  {
    return _begin[range] + indexInRange(tid, range);
  }
  /**
   * Get the position of a thread within its range
   * @param tid The thread
   * @param range The range of the thread, from range()
   */
  KOKKOS_FUNCTION MOOSE_KOKKOS_INDEX_TYPE indexInRange(const MOOSE_KOKKOS_INDEX_TYPE tid,
                                                       const unsigned int range) const
  {
    return tid - _offset[range];
  }
#endif

private:
  /// Number of faces
  MOOSE_KOKKOS_INDEX_TYPE _size = 0;
  /// Number of ranges
  unsigned int _num_ranges = 0;
  /// First face of each range
  Array<ContiguousFaceID> _begin;
  /// First thread of each range, and the number of faces at the end
  Array<MOOSE_KOKKOS_INDEX_TYPE> _offset;
};

/**
 * The Kokkos interface that holds the host reference of the linear finite volume mesh data and
 * copies it to device during parallel dispatch
 */
class FVMeshHolder
{
public:
  /**
   * Constructor
   * @param mesh The linear finite volume mesh data
   */
  FVMeshHolder(const FVMesh & mesh) : _fv_mesh_host(mesh), _fv_mesh_device(mesh) {}
  /**
   * Copy constructor
   */
  FVMeshHolder(const FVMeshHolder & holder)
    : _fv_mesh_host(holder._fv_mesh_host), _fv_mesh_device(holder._fv_mesh_host)
  {
  }

#ifdef MOOSE_KOKKOS_SCOPE
  /**
   * Get the linear finite volume mesh data on the architecture this function is called on
   */
  KOKKOS_FUNCTION const FVMesh & kokkosFVMesh() const
  {
    KOKKOS_IF_ON_HOST(
        if (!_fv_mesh_host.initialized()) mooseError(
            "kokkosFVMesh() was called too early. The Kokkos finite volume mesh data is available "
            "after problem initialization. Override initialSetup() if you need to setup your "
            "object data using it.");

        return _fv_mesh_host;)

    return _fv_mesh_device;
  }
#endif

private:
  /// Host reference of the linear finite volume mesh data
  const FVMesh & _fv_mesh_host;
  /// Device copy of the linear finite volume mesh data
  const FVMesh _fv_mesh_device;
};

} // namespace Moose::Kokkos
