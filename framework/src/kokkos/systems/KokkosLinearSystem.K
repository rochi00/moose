//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearSystem.h"

#include "Attributes.h"
#include "KokkosLinearFVElementalKernel.h"
#include "KokkosLinearFVFluxKernel.h"
#include "KokkosLinearFVBoundaryCondition.h"
#include "KokkosFVSystemAssembly.h"
#include "MooseUtils.h"
#include "MaterialPropertyInterface.h"
#include "MaterialBase.h"
#include "FEProblemBase.h"

#include "libmesh/equation_systems.h"

#include <algorithm>
#include <optional>
#include <set>
#include <unordered_set>

namespace
{
/// The linear finite volume kernels and boundary conditions of a system that are, or are not,
/// Kokkos objects
std::vector<MooseObject *>
linearFVObjects(FEProblemBase & problem, const unsigned int sys_num, const bool kokkos)
{
  std::vector<MooseObject *> all_objects;
  for (const auto & system :
       {"LinearFVElementalKernel", "LinearFVFluxKernel", "LinearFVBoundaryCondition"})
  {
    std::vector<MooseObject *> objects;
    problem.theWarehouse()
        .query()
        .template condition<AttribSysNum>(sys_num)
        .template condition<AttribSystem>(system)
        .template condition<AttribKokkos>(kokkos)
        .template condition<AttribThread>(0)
        .queryInto(objects);
    all_objects.insert(all_objects.end(), objects.begin(), objects.end());
  }
  return all_objects;
}
}

void
LinearSystem::addKokkosKernel(const std::string & kernel_name,
                              const std::string & name,
                              InputParameters & parameters)
{
  auto kernel = _factory.create<Moose::Kokkos::LinearFVKernel>(kernel_name, name, parameters);

  if (!std::dynamic_pointer_cast<Moose::Kokkos::LinearFVElementalKernel>(kernel) &&
      !std::dynamic_pointer_cast<Moose::Kokkos::LinearFVFluxKernel>(kernel))
    mooseError("Unknown Kokkos LinearFVKernel type for object named ", kernel->name());

  _fe_problem.theWarehouse().add(kernel);
}

void
LinearSystem::addKokkosBoundaryCondition(const std::string & bc_name,
                                         const std::string & name,
                                         InputParameters & parameters)
{
  auto bc = _factory.create<Moose::Kokkos::LinearFVBoundaryCondition>(bc_name, name, parameters);
  _fe_problem.theWarehouse().add(bc);

  const auto & boundary_ids = bc->boundaryIDs();
  const auto * const bc_var = bc->variable().mooseFieldVar();
  _vars[0].addBoundaryVar(boundary_ids, bc_var);
}

void
LinearSystem::selectKokkosPetscTypes()
{
  if (linearFVObjects(_fe_problem, number(), true).empty() ||
      !linearFVObjects(_fe_problem, number(), false).empty())
    return;

  // Added as defaults to the problem PETSc options, which are applied to the existing vectors and
  // matrices when the options are inserted before the first solve. A Kokkos matrix needs Kokkos
  // vectors, so nothing is added if the input file (stored upper case) or the command line sets
  // this system's matrix type or another vector type.
  auto & options = _fe_problem.getPetscOptions().pairs;
  auto user_value = [this, &options](const std::string & option) -> std::optional<std::string>
  {
    char value[256] = "";
    PetscBool on_command_line = PETSC_FALSE;
    LibmeshPetscCall(PetscOptionsGetString(
        nullptr, nullptr, option.c_str(), value, sizeof(value), &on_command_line));
    if (on_command_line)
      return std::string(value);

    for (const auto & [name, pair_value] : options)
      if (MooseUtils::toLower(name) == option)
        return pair_value;

    return {};
  };

  // MOOSE applies matrix types only through options prefixed by the solver system name
  const auto mat_type_option = "-" + name() + "_mat_type";
  const auto vec_type = user_value("-vec_type");
  if (user_value(mat_type_option) || (vec_type && *vec_type != "kokkos"))
    return;

  options.emplace_back(mat_type_option, "aijkokkos");
  // Another Kokkos linear system may have added it already
  if (!vec_type)
    options.emplace_back("-vec_type", "kokkos");

  _console << "Linear system '" << name()
           << "' is assembled by Kokkos objects only; keeping its matrix and vectors in device "
              "memory with the PETSc options: "
           << mat_type_option << " aijkokkos -vec_type kokkos" << std::endl;
}

void
LinearSystem::initialSetupKokkosLinearFV()
{
  // Kokkos objects are currently instantiated only for thread 0 (single copy per object).
  // All three queries below use AttribThread(0) for this reason. If Kokkos objects are ever
  // allowed to be multi-threaded, this function will need to loop over threads and each kernel's
  // initialSetup will work correctly via its own _tid member.
  std::vector<Moose::Kokkos::LinearFVElementalKernel *> elemental_kernels;
  _fe_problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(number())
      .template condition<AttribSystem>("LinearFVElementalKernel")
      .template condition<AttribKokkos>(true)
      .template condition<AttribThread>(0)
      .queryInto(elemental_kernels);

  for (auto * kernel : elemental_kernels)
    kernel->initialSetup();

  std::vector<Moose::Kokkos::LinearFVBoundaryCondition *> bcs;
  _fe_problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(number())
      .template condition<AttribSystem>("LinearFVBoundaryCondition")
      .template condition<AttribKokkos>(true)
      .template condition<AttribThread>(0)
      .queryInto(bcs);

  for (auto * bc : bcs)
    bc->initialSetup();

  std::vector<Moose::Kokkos::LinearFVFluxKernel *> flux_kernels;
  _fe_problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(number())
      .template condition<AttribSystem>("LinearFVFluxKernel")
      .template condition<AttribKokkos>(true)
      .template condition<AttribThread>(0)
      .queryInto(flux_kernels);

  for (auto * kernel : flux_kernels)
    kernel->initialSetup();

  // A system with Kokkos linear FV objects must not also contain host linear FV objects: those
  // force the host assembly loops over every element and face to run in addition to the Kokkos
  // assembly. Unless explicitly allowed, refuse such systems instead of silently running slower.
  const auto host_objects = linearFVObjects(_fe_problem, number(), false);

  const bool has_kokkos_objects =
      !elemental_kernels.empty() || !flux_kernels.empty() || !bcs.empty();

  if (has_kokkos_objects && !host_objects.empty() &&
      !_fe_problem.getParam<bool>("allow_mixed_kokkos_linear_fv"))
  {
    std::string list;
    for (const auto * object : host_objects)
      list += "\n  '" + object->name() + "' (type '" + object->type() + "')";
    mooseError("Linear system '",
               name(),
               "' contains Kokkos linear finite volume objects together with the following "
               "non-Kokkos linear finite volume objects:",
               list,
               "\nMixing requires the host assembly loops to run alongside the Kokkos assembly at "
               "every evaluation. Replace them with Kokkos objects, or set "
               "'Problem/allow_mixed_kokkos_linear_fv = true' to accept the slower mixed "
               "assembly.");
  }

  _all_linear_fv_objects_kokkos = has_kokkos_objects && host_objects.empty();
}

void
LinearSystem::computeKokkosLinearSystem(const std::set<TagID> & vector_tags,
                                        const std::set<TagID> & matrix_tags)
{
  TIME_SECTION("computeKokkosLinearSystem", 1);

  std::vector<Moose::Kokkos::LinearFVElementalKernel *> elemental_for_vector_tags;
  _fe_problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(number())
      .template condition<AttribSystem>("LinearFVElementalKernel")
      .template condition<AttribKokkos>(true)
      .template condition<AttribThread>(0)
      .template condition<AttribVectorTags>(vector_tags)
      .queryInto(elemental_for_vector_tags);

  std::vector<Moose::Kokkos::LinearFVElementalKernel *> elemental_for_matrix_tags;
  _fe_problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(number())
      .template condition<AttribSystem>("LinearFVElementalKernel")
      .template condition<AttribKokkos>(true)
      .template condition<AttribThread>(0)
      .template condition<AttribMatrixTags>(matrix_tags)
      .queryInto(elemental_for_matrix_tags);

  std::vector<Moose::Kokkos::LinearFVFluxKernel *> flux_for_vector_tags;
  _fe_problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(number())
      .template condition<AttribSystem>("LinearFVFluxKernel")
      .template condition<AttribKokkos>(true)
      .template condition<AttribThread>(0)
      .template condition<AttribVectorTags>(vector_tags)
      .queryInto(flux_for_vector_tags);

  std::vector<Moose::Kokkos::LinearFVFluxKernel *> flux_for_matrix_tags;
  _fe_problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(number())
      .template condition<AttribSystem>("LinearFVFluxKernel")
      .template condition<AttribKokkos>(true)
      .template condition<AttribThread>(0)
      .template condition<AttribMatrixTags>(matrix_tags)
      .queryInto(flux_for_matrix_tags);

  std::vector<Moose::Kokkos::LinearFVBoundaryCondition *> bcs;
  _fe_problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(number())
      .template condition<AttribSystem>("LinearFVBoundaryCondition")
      .template condition<AttribKokkos>(true)
      .template condition<AttribThread>(0)
      .queryInto(bcs);

  std::vector<Moose::Kokkos::LinearFVElementalKernel *> elemental_kernels;
  MooseUtils::getUnion(elemental_for_vector_tags, elemental_for_matrix_tags, elemental_kernels);

  std::vector<Moose::Kokkos::LinearFVFluxKernel *> flux_kernels;
  MooseUtils::getUnion(flux_for_vector_tags, flux_for_matrix_tags, flux_kernels);

  if (elemental_kernels.empty() && flux_kernels.empty())
    return;

  auto & our_kokkos_system = _fe_problem.getKokkosSystem(number());

  std::set<TagID> needed_solution_vector_tags;

  for (const auto & tag : _fe_problem.getVectorTags(Moose::VECTOR_TAG_SOLUTION))
    needed_solution_vector_tags.insert(tag._id);

  // Kokkos linear FV kernels and boundary conditions read only their own variable, which belongs to
  // this system, so the other Kokkos systems are not synchronized. The Kokkos system provides the
  // solution vectors; the tagged matrices and vectors being assembled bypass it.
  auto & kokkos_systems = _fe_problem.getKokkosSystems();
  if (!_kokkos_fv_assembly)
    _kokkos_fv_assembly = std::make_unique<Moose::Kokkos::FVSystemAssembly>(*this, our_kokkos_system);
  {
    TIME_SECTION("KokkosFVSyncToDevice", 2);
    our_kokkos_system.setActiveSolutionTags(needed_solution_vector_tags);
    our_kokkos_system.sync(Moose::Kokkos::MemcpyType::HOST_TO_DEVICE);

    kokkos_systems.copyToDevice();

    _kokkos_fv_assembly->begin(vector_tags, matrix_tags);
  }

  std::unordered_set<unsigned int> value_bc_vars;
  std::unordered_set<unsigned int> normal_gradient_bc_vars;

  for (const auto kernel : flux_kernels)
  {
    if (kernel->needsBoundaryValueData())
      value_bc_vars.insert(kernel->variable().var());
    if (kernel->needsBoundaryNormalGradientData())
      normal_gradient_bc_vars.insert(kernel->variable().var());
  }

  // The boundary data and assembly loops are dispatched without fences. They run in order on the
  // default execution space instance, so each loop sees the results of the loops before it. On a
  // host backend the dispatches are synchronous and the sections below measure execution; on a
  // device they measure launch, and KokkosFVWait the remaining execution.
  {
    TIME_SECTION("KokkosFVBoundaryData", 2);
    for (auto bc : bcs)
    {
      const auto bc_var = bc->variable().var();
      if (value_bc_vars.count(bc_var))
        bc->computeBoundaryValueData();
      if (normal_gradient_bc_vars.count(bc_var))
        bc->computeBoundaryNormalGradientData();
    }
  }

  // Each kernel adds its matrix and right-hand side contributions in one pass: to the matrices if
  // its matrix tags are requested, and to the vectors if its vector tags are
  auto contains = [](const auto & kernels, const auto * const kernel)
  { return std::find(kernels.begin(), kernels.end(), kernel) != kernels.end(); };

  {
    TIME_SECTION("KokkosFVElemental", 2);
    for (auto kernel : elemental_kernels)
      kernel->computeContributions(contains(elemental_for_matrix_tags, kernel),
                                   contains(elemental_for_vector_tags, kernel),
                                   *_kokkos_fv_assembly);
  }
  {
    TIME_SECTION("KokkosFVFlux", 2);
    for (auto kernel : flux_kernels)
      kernel->computeContributions(contains(flux_for_matrix_tags, kernel),
                                   contains(flux_for_vector_tags, kernel),
                                   *_kokkos_fv_assembly);
  }

  {
    TIME_SECTION("KokkosFVWait", 2);
    // The only synchronization of the assembly: the results are read below
    ::Kokkos::fence();
  }

  {
    TIME_SECTION("KokkosFVSyncToHost", 2);
    our_kokkos_system.sync(Moose::Kokkos::MemcpyType::DEVICE_TO_HOST);
  }
  {
    TIME_SECTION("KokkosFVUpload", 2);
    _kokkos_fv_assembly->end();
  }

  our_kokkos_system.clearActiveVariables();
  our_kokkos_system.clearActiveSolutionTags();
}
