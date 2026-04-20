//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "libmesh/gpu/kokkos_fe_types.h"

#include "libmesh/libmesh_common.h"

namespace libMesh::Kokkos
{

FEElemTopology
toKokkosTopology(libMesh::ElemType t)
{
  switch (t)
  {
    case libMesh::EDGE2:  return FEElemTopology::EDGE2;
    case libMesh::EDGE3:  return FEElemTopology::EDGE3;
    case libMesh::TRI3:   return FEElemTopology::TRI3;
    case libMesh::TRI6:   return FEElemTopology::TRI6;
    case libMesh::QUAD4:  return FEElemTopology::QUAD4;
    case libMesh::QUAD8:  return FEElemTopology::QUAD8;
    case libMesh::QUAD9:  return FEElemTopology::QUAD9;
    case libMesh::TET4:   return FEElemTopology::TET4;
    case libMesh::TET10:  return FEElemTopology::TET10;
    case libMesh::HEX8:   return FEElemTopology::HEX8;
    case libMesh::HEX20:  return FEElemTopology::HEX20;
    case libMesh::HEX27:  return FEElemTopology::HEX27;
    default:
      libmesh_error_msg("Element type " << t << " is not supported by the Kokkos native FE path");
  }
}

FEFamily
toKokkosFamily(libMesh::FEFamily f)
{
  switch (f)
  {
    case libMesh::LAGRANGE:     return FEFamily::LAGRANGE;
    case libMesh::LAGRANGE_VEC: return FEFamily::LAGRANGE_VEC;
    case libMesh::HERMITE:      return FEFamily::HERMITE;
    case libMesh::MONOMIAL:     return FEFamily::MONOMIAL;
    case libMesh::MONOMIAL_VEC: return FEFamily::MONOMIAL_VEC;
    default:
      // Return a sentinel rather than erroring — the shape-key table is populated
      // for every FE type in the problem, so an unrecognised family must not abort.
      // The family_supported gate in initShape() will route these types to the
      // libMesh phi-table fallback path.
      return FEFamily::UNKNOWN;
  }
}

FEElemTopology
getSideTopology(FEElemTopology parent)
{
  switch (parent)
  {
    // 1D: sides are vertex nodes, represented as EDGE2 (degenerate)
    case FEElemTopology::EDGE2:
    case FEElemTopology::EDGE3:
      return FEElemTopology::EDGE2;

    // 2D first-order: sides are linear edges
    case FEElemTopology::TRI3:
    case FEElemTopology::QUAD4:
      return FEElemTopology::EDGE2;

    // 2D second-order: sides are quadratic edges
    case FEElemTopology::TRI6:
    case FEElemTopology::QUAD8:
    case FEElemTopology::QUAD9:
      return FEElemTopology::EDGE3;

    // 3D first-order: sides are linear triangles / quads
    case FEElemTopology::TET4:
      return FEElemTopology::TRI3;
    case FEElemTopology::HEX8:
      return FEElemTopology::QUAD4;

    // 3D second-order: sides are quadratic triangles / quads
    case FEElemTopology::TET10:
      return FEElemTopology::TRI6;
    case FEElemTopology::HEX20:
      return FEElemTopology::QUAD8;
    case FEElemTopology::HEX27:
      return FEElemTopology::QUAD9;

    default:
      libmesh_error_msg("getSideTopology: unknown FEElemTopology value "
                 << static_cast<unsigned int>(parent));
  }
}

unsigned int
nDofs(FEFamily family, FEElemTopology topo)
{
  switch (family)
  {
    case FEFamily::LAGRANGE:
    {
      switch (topo)
      {
        case FEElemTopology::EDGE2:  return 2;
        case FEElemTopology::EDGE3:  return 3;
        case FEElemTopology::TRI3:   return 3;
        case FEElemTopology::TRI6:   return 6;
        case FEElemTopology::QUAD4:  return 4;
        case FEElemTopology::QUAD8:  return 8;
        case FEElemTopology::QUAD9:  return 9;
        case FEElemTopology::TET4:   return 4;
        case FEElemTopology::TET10:  return 10;
        case FEElemTopology::HEX8:   return 8;
        case FEElemTopology::HEX20:  return 20;
        case FEElemTopology::HEX27:  return 27;
        default:
          libmesh_error_msg("nDofs: unsupported topology " << static_cast<unsigned int>(topo)
                     << " for LAGRANGE family");
      }
    }

    case FEFamily::HERMITE:
      libmesh_error_msg("nDofs: HERMITE family is not supported by the native FE path; "
                 "use the libMesh fallback path instead");

    case FEFamily::LAGRANGE_VEC:
      libmesh_error_msg("nDofs: LAGRANGE_VEC requires a spatial dimension argument; "
                 "not yet supported by the native FE path");

    default:
      libmesh_error_msg("nDofs: unsupported FE family " << static_cast<unsigned int>(family));
  }
}

FEElemClass
classFromTopology(FEElemTopology topo)
{
  switch (topo)
  {
    case FEElemTopology::EDGE2:
    case FEElemTopology::EDGE3:
      return FEElemClass::EDGE;

    case FEElemTopology::TRI3:
    case FEElemTopology::TRI6:
      return FEElemClass::TRI;

    case FEElemTopology::QUAD4:
    case FEElemTopology::QUAD8:
    case FEElemTopology::QUAD9:
      return FEElemClass::QUAD;

    case FEElemTopology::TET4:
    case FEElemTopology::TET10:
      return FEElemClass::TET;

    case FEElemTopology::HEX8:
    case FEElemTopology::HEX20:
    case FEElemTopology::HEX27:
      return FEElemClass::HEX;

    default:
      libmesh_error_msg("classFromTopology: unsupported FEElemTopology value "
                 << static_cast<unsigned int>(topo));
  }
}

} // namespace libMesh::Kokkos
