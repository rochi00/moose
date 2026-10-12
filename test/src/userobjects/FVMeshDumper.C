//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVMeshDumper.h"
#include "MooseMesh.h"
#include "MooseVariableFieldBase.h"
#include "SystemBase.h"
#include "FaceInfo.h"
#include "ElemInfo.h"

#include <fstream>

registerMooseObject("MooseTestApp", FVMeshDumper);

InputParameters
FVMeshDumper::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription("Writes finite volume cell and face data to a binary file for the "
                             "standalone Kokkos linear FV layout benchmark.");
  params.addRequiredParam<FileName>("file", "Output file name");
  params.addRequiredParam<VariableName>("variable",
                                        "Linear FV variable whose dof indices are written");
  params.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;
  return params;
}

FVMeshDumper::FVMeshDumper(const InputParameters & params) : GeneralUserObject(params) {}

namespace
{
template <typename T>
void
put(std::ofstream & out, const T & value)
{
  out.write(reinterpret_cast<const char *>(&value), sizeof(T));
}

void
putPoint(std::ofstream & out, const Point & p)
{
  for (const auto d : make_range(3))
    put<double>(out, p(d));
}
}

void
FVMeshDumper::execute()
{
  if (n_processors() > 1)
    mooseError("FVMeshDumper only supports serial runs");

  const auto & var = _fe_problem.getVariable(0, getParam<VariableName>("variable"));
  const auto sys_num = var.sys().number();
  const auto var_num = var.number();

  // Cells in libMesh element id order, which is also the native Kokkos contiguous order for a
  // single-subdomain serial mesh
  const auto & mesh = _fe_problem.mesh();
  std::vector<const ElemInfo *> cells(mesh.elemInfoVector().begin(), mesh.elemInfoVector().end());
  std::sort(cells.begin(),
            cells.end(),
            [](const ElemInfo * a, const ElemInfo * b)
            { return a->elem()->id() < b->elem()->id(); });
  std::unordered_map<dof_id_type, int64_t> cell_index;
  for (const auto i : index_range(cells))
    cell_index[cells[i]->elem()->id()] = i;

  const auto & faces = mesh.faceInfo();

  std::ofstream out(getParam<FileName>("file"), std::ios::binary);
  if (!out)
    mooseError("Cannot open ", getParam<FileName>("file"));

  // Header: magic, version, counts
  out.write("FVD1", 4);
  put<int64_t>(out, cells.size());
  put<int64_t>(out, faces.size());

  // Cell records: libMesh element id, dof index, volume (with coordinate factor), centroid
  for (const auto * ei : cells)
  {
    put<int64_t>(out, ei->elem()->id());
    put<int64_t>(out, ei->dofIndices()[sys_num][var_num]);
    put<double>(out, ei->volume() * ei->coordFactor());
    putPoint(out, ei->centroid());
  }

  // Face records: elem and neighbor cell indices (-1 if none), sides, area (with coordinate
  // factor), normal, face centroid, dCN, gC, boundary id (-1 if none)
  for (const auto * fi : faces)
  {
    put<int64_t>(out, libmesh_map_find(cell_index, fi->elem().id()));
    put<int64_t>(out, fi->neighborPtr() ? libmesh_map_find(cell_index, fi->neighbor().id()) : -1);
    put<int32_t>(out, fi->elemSideID());
    put<int32_t>(out, fi->neighborPtr() ? fi->neighborSideID() : -1);
    put<double>(out, fi->faceArea() * fi->faceCoord());
    putPoint(out, fi->normal());
    putPoint(out, fi->faceCentroid());
    putPoint(out, fi->dCN());
    put<double>(out, fi->gC());
    put<int32_t>(out, fi->boundaryIDs().empty() ? -1 : int32_t(*fi->boundaryIDs().begin()));
  }
}
