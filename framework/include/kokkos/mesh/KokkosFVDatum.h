//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosFVMesh.h"

namespace Moose::Kokkos
{

/**
 * The face seen by a Kokkos linear finite volume object, built per thread. The elem cell is the
 * cell the face is evaluated for and the normal points out of it. It is the elem cell of the
 * FaceInfo unless the datum is flipped, which a boundary face of a variable that lives only on
 * the neighbor side of the FaceInfo needs.
 */
class FaceDatum
{
public:
  /**
   * Constructor
   * @param face The contiguous face ID
   * @param flip Whether the face is seen from the neighbor cell of its FaceInfo
   * @param mesh The Kokkos finite volume mesh data
   */
  KOKKOS_FUNCTION FaceDatum(const ContiguousFaceID face, const bool flip, const FVMesh & mesh)
    : _face(face), _flip(flip), _mesh(mesh)
  {
  }

  /**
   * Get the contiguous face ID
   */
  KOKKOS_FUNCTION ContiguousFaceID face() const { return _face; }
  /**
   * Get the contiguous ID of the cell this face is evaluated for
   */
  KOKKOS_FUNCTION ContiguousElementID elem() const
  {
    return _flip ? _mesh.getNeighbor(_face) : _mesh.getElem(_face);
  }
  /**
   * Get the contiguous ID of the cell on the other side of the face
   * @returns The contiguous ID, or libMesh::DofObject::invalid_id on a mesh boundary
   */
  KOKKOS_FUNCTION ContiguousElementID neighbor() const
  {
    return _flip ? _mesh.getElem(_face) : _mesh.getNeighbor(_face);
  }
  /**
   * Get the coordinate-weighted face area
   */
  KOKKOS_FUNCTION Real faceArea() const { return _mesh.getArea(_face); }
  /**
   * Get the unit face normal, pointing out of elem()
   */
  KOKKOS_FUNCTION Real3 faceNormal() const
  {
    return _flip ? -_mesh.getNormal(_face) : _mesh.getNormal(_face);
  }
  /**
   * Get the face centroid
   */
  KOKKOS_FUNCTION Real3 faceCentroid() const { return _mesh.getCentroid(_face); }
  /**
   * Get the distance between the centroids of the two cells, or between the cell centroid and the
   * face centroid on a mesh boundary
   */
  KOKKOS_FUNCTION Real faceDCNMag() const { return _mesh.getDCNMag(_face); }
  /**
   * Get the distance between the centroid of elem() and the face centroid
   */
  KOKKOS_FUNCTION Real faceDCFMag() const
  {
    return (_mesh.getCentroid(_face) - _mesh.getElemCentroid(elem())).norm();
  }
  /**
   * Get the centroid of elem()
   */
  KOKKOS_FUNCTION Real3 elemCentroid() const { return _mesh.getElemCentroid(elem()); }
  /**
   * Get the coordinate-weighted volume of elem()
   */
  KOKKOS_FUNCTION Real elemVolume() const { return _mesh.getElemVolume(elem()); }

private:
  /// Contiguous face ID
  const ContiguousFaceID _face;
  /// Whether the face is seen from the neighbor cell of its FaceInfo
  const bool _flip;
  /// Kokkos finite volume mesh data
  const FVMesh & _mesh;
};

/**
 * The cell seen by a Kokkos linear finite volume elemental kernel, built per thread
 */
class CellDatum
{
public:
  /**
   * Constructor
   * @param elem The contiguous ID of an owned cell
   * @param mesh The Kokkos finite volume mesh data
   */
  KOKKOS_FUNCTION CellDatum(const ContiguousElementID elem, const FVMesh & mesh)
    : _elem(elem), _mesh(mesh)
  {
  }

  /**
   * Get the contiguous cell ID
   */
  KOKKOS_FUNCTION ContiguousElementID elem() const { return _elem; }
  /**
   * Get the coordinate-weighted cell volume
   */
  KOKKOS_FUNCTION Real elemVolume() const { return _mesh.getElemVolume(_elem); }
  /**
   * Get the cell centroid
   */
  KOKKOS_FUNCTION Real3 elemCentroid() const { return _mesh.getElemCentroid(_elem); }

private:
  /// Contiguous cell ID
  const ContiguousElementID _elem;
  /// Kokkos finite volume mesh data
  const FVMesh & _mesh;
};

} // namespace Moose::Kokkos

using FaceDatum = Moose::Kokkos::FaceDatum;
using CellDatum = Moose::Kokkos::CellDatum;
