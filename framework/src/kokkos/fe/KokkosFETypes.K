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

unsigned int
nDofs(FEFamily family, libMesh::ElemType topo)
{
  switch (family)
  {
    case FEFamily::LAGRANGE:
    {
      switch (topo)
      {
        case libMesh::EDGE2:  return 2;
        case libMesh::EDGE3:  return 3;
        case libMesh::TRI3:   return 3;
        case libMesh::TRI6:   return 6;
        case libMesh::QUAD4:  return 4;
        case libMesh::QUAD8:  return 8;
        case libMesh::QUAD9:  return 9;
        case libMesh::TET4:   return 4;
        case libMesh::TET10:  return 10;
        case libMesh::HEX8:   return 8;
        case libMesh::HEX20:  return 20;
        case libMesh::HEX27:  return 27;
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

} // namespace libMesh::Kokkos
