//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConservativeDiffuseInterfaceRhieChowMassFlux.h"

registerMooseObject("NavierStokesApp", ConservativeDiffuseInterfaceRhieChowMassFlux);
#include "NSFVLeastSquares.h"
#include "DerivedFaceFunctor.h"

#include "LinearFVAdvectionDiffusionBC.h"
#include "LinearFVFluxKernel.h"
#include "LinearFVPressureInletOutletVelocityBC.h"
#include "LinearSystem.h"
#include "MooseFunctorArguments.h"
#include "MooseMesh.h"
#include "SubProblem.h"
#include "FVUtils.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>

namespace
{
constexpr Real pressure_writeback_face_ainv_relative_tolerance = 1e-4;

template <typename FaceNormalScalar>
RealVectorValue
reconstructCellVectorFromFaceNormalScalars(const MooseMesh & moose_mesh,
                                           const unsigned int dim,
                                           const ElemInfo * elem_info,
                                           FaceNormalScalar face_normal_scalar)
{
  if (!elem_info)
    return RealVectorValue();

  const Elem * const elem = elem_info->elem();
  if (!elem)
    return RealVectorValue();

  DenseMatrix<Real> normal_matrix(dim, dim);
  DenseVector<Real> rhs(dim);
  DenseVector<Real> solution(dim);

  for (const auto i : make_range(dim))
  {
    rhs(i) = 0.0;
    solution(i) = 0.0;
    for (const auto j : make_range(dim))
      normal_matrix(i, j) = 0.0;
  }

  for (const auto side : make_range(elem->n_sides()))
  {
    const FaceInfo * const fi_loc = NS::faceInfoForElementSide(moose_mesh, elem, side);
    if (!fi_loc)
      continue;

    const Real face_measure = NS::faceMeasure(*fi_loc);
    if (face_measure <= 0.0)
      continue;

    const Real orientation = fi_loc->elemPtr() == elem ? 1.0 : -1.0;
    Real normal_scalar = 0.0;
    if (!face_normal_scalar(*fi_loc, orientation, normal_scalar))
      continue;

    const RealVectorValue outward_normal = orientation * fi_loc->normal();
    const Real psi_f = orientation * normal_scalar;
    for (const auto i : make_range(dim))
    {
      rhs(i) += face_measure * outward_normal(i) * psi_f;
      for (const auto j : make_range(dim))
        normal_matrix(i, j) += face_measure * outward_normal(i) * outward_normal(j);
    }
  }

  bool singular = false;
  if (dim == 1)
    singular = std::abs(normal_matrix(0, 0)) <= libMesh::TOLERANCE;
  else if (dim == 2)
    singular = std::abs(normal_matrix(0, 0) * normal_matrix(1, 1) -
                        normal_matrix(0, 1) * normal_matrix(1, 0)) <= libMesh::TOLERANCE;
  else if (dim == 3)
  {
    const Real det = normal_matrix(0, 0) * (normal_matrix(1, 1) * normal_matrix(2, 2) -
                                            normal_matrix(1, 2) * normal_matrix(2, 1)) -
                     normal_matrix(0, 1) * (normal_matrix(1, 0) * normal_matrix(2, 2) -
                                            normal_matrix(1, 2) * normal_matrix(2, 0)) +
                     normal_matrix(0, 2) * (normal_matrix(1, 0) * normal_matrix(2, 1) -
                                            normal_matrix(1, 1) * normal_matrix(2, 0));
    singular = std::abs(det) <= libMesh::TOLERANCE;
  }

  if (singular)
    return RealVectorValue();

  normal_matrix.lu_solve(rhs, solution);

  RealVectorValue reconstructed_source;
  for (const auto i : make_range(dim))
    reconstructed_source(i) = solution(i);

  return reconstructed_source;
}
}

struct ConservativeDiffuseInterfaceRhieChowMassFlux::MomentumPredictorOperator
{
  Real relaxation_parameter = 1.0;
  std::unique_ptr<NumericVector<Number>> rhs_raw;
  std::unique_ptr<NumericVector<Number>> pressure_coupling_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> h_diagonal_correction_raw;
  std::unique_ptr<NumericVector<Number>> interior_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> boundary_component_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> boundary_relax_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> boundary_a_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> boundary_dominance_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> offdiag_abs_sum_raw;
  struct BoundaryMatrixContribution
  {
    dof_id_type face_id;
    dof_id_type dof;
    Real contribution;
  };
  std::vector<BoundaryMatrixContribution> boundary_matrix_contributions;
};

const ConservativeDiffuseInterfaceRhieChowMassFlux::MomentumPredictorOperator *
ConservativeDiffuseInterfaceRhieChowMassFlux::cachedMomentumPredictorOperator(
    const unsigned int system_i) const
{
  if (system_i >= _momentum_predictor_operators.size())
    return nullptr;

  return _momentum_predictor_operators[system_i].get();
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::beginMomentumPredictorAssembly(
    const unsigned int system_i)
{
  // The first component starts a new vector operator.
  if (system_i == 0)
  {
    _momentum_predictor_operators.clear();
    _predictor_ready = false;
  }

  mooseAssert(system_i < _momentum_implicit_systems.size() && _momentum_implicit_systems[system_i],
              "The requested momentum component is not linked to RhieChowMassFlux.");

  if (_momentum_predictor_operators.size() != _momentum_systems.size())
    _momentum_predictor_operators.resize(_momentum_systems.size());

  auto * momentum_system = _momentum_implicit_systems[system_i];
  NumericVector<Number> & current_local_solution = *(momentum_system->current_local_solution);
  auto predictor = std::make_unique<MomentumPredictorOperator>();

  auto initialize_zero_vector =
      [&current_local_solution](std::unique_ptr<NumericVector<Number>> & v)
  {
    v = current_local_solution.zero_clone();
    v->zero();
    v->close();
  };

  for (auto * vector : {&predictor->rhs_raw,
                        &predictor->pressure_coupling_diagonal_raw,
                        &predictor->h_diagonal_correction_raw,
                        &predictor->interior_diagonal_raw,
                        &predictor->boundary_component_diagonal_raw,
                        &predictor->boundary_relax_diagonal_raw,
                        &predictor->boundary_a_diagonal_raw,
                        &predictor->boundary_dominance_diagonal_raw,
                        &predictor->offdiag_abs_sum_raw})
    initialize_zero_vector(*vector);

  _momentum_predictor_operators[system_i] = std::move(predictor);
  _predictor_ready = false;
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::momentumPredictorAssembled(
    const unsigned int system_i,
    SparseMatrix<Number> & matrix,
    const NumericVector<Number> & rhs,
    const Real relaxation_parameter)
{
  mooseAssert(system_i < _momentum_predictor_operators.size() &&
                  _momentum_predictor_operators[system_i],
              "A split momentum predictor operator assembly was not started.");

  auto * mmat = dynamic_cast<PetscMatrix<Number> *>(&matrix);
  mooseAssert(mmat, "The split momentum predictor cache requires a PETSc-backed momentum matrix.");

  auto & predictor = *_momentum_predictor_operators[system_i];

  *predictor.rhs_raw = rhs;
  predictor.rhs_raw->close();

  // The assembled diagonal mixes interior coefficients with the boundary ones, which fvMatrix
  // treats differently (component-averaged for A, per component for H). They are not recoverable
  // from the matrix, so ask the flux kernels that assembled them, exactly as the pressure
  // diffusion kernel is asked for the pressure face flux.
  predictor.boundary_matrix_contributions.clear();
  std::vector<LinearFVFluxKernel *> flux_kernels;
  _fe_problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(_global_momentum_system_numbers[system_i])
      .template condition<AttribSystem>("LinearFVFluxKernel")
      .template condition<AttribThread>(_tid)
      .queryInto(flux_kernels);

  for (const auto * fi : _flow_face_info)
    for (auto * kernel : flux_kernels)
    {
      const bool on_elem = kernel->hasBlocks(fi->elemSubdomainID());
      const bool on_neighbor = fi->neighborPtr() && kernel->hasBlocks(fi->neighborSubdomainID());
      if (!on_elem && !on_neighbor)
        continue;

      auto & var = kernel->variable();
      const auto face_type = fi->faceType(std::make_pair(var.number(), var.sys().number()));
      if (face_type != FaceInfo::VarFaceNeighbors::ELEM &&
          face_type != FaceInfo::VarFaceNeighbors::NEIGHBOR)
        continue;

      auto * const bc = var.getBoundaryCondition(*fi);
      if (!bc)
        continue;

      kernel->setupFaceData(fi);
      kernel->setCurrentFaceArea(NS::faceMeasure(*fi));
      const Real contribution = kernel->computeBoundaryMatrixContribution(*bc);
      const auto & elem_info =
          face_type == FaceInfo::VarFaceNeighbors::ELEM ? *fi->elemInfo() : *fi->neighborInfo();
      predictor.boundary_matrix_contributions.push_back(
          {fi->id(), elem_info.dofIndices()[var.sys().number()][var.number()], contribution});
    }

  predictor.interior_diagonal_raw->zero();
  predictor.interior_diagonal_raw->close();
  mmat->get_diagonal(*predictor.interior_diagonal_raw);
  predictor.interior_diagonal_raw->close();
  for (const auto & contribution : predictor.boundary_matrix_contributions)
    predictor.interior_diagonal_raw->add(contribution.dof, -contribution.contribution);
  predictor.interior_diagonal_raw->close();

  predictor.offdiag_abs_sum_raw->zero();
  predictor.offdiag_abs_sum_raw->close();

  std::vector<dof_id_type> offdiag_abs_sum_indices;
  std::vector<Real> offdiag_abs_sum_values;
  offdiag_abs_sum_indices.reserve(matrix.row_stop() - matrix.row_start());
  offdiag_abs_sum_values.reserve(matrix.row_stop() - matrix.row_start());

  for (const auto row_dof : make_range(matrix.row_start(), matrix.row_stop()))
  {
    std::vector<numeric_index_type> column_dofs;
    std::vector<Real> coefficients;
    mmat->get_row(row_dof, column_dofs, coefficients);

    Real offdiag_abs_sum = 0.0;
    for (const auto entry_i : index_range(column_dofs))
      if (column_dofs[entry_i] != row_dof)
        offdiag_abs_sum += std::abs(coefficients[entry_i]);

    offdiag_abs_sum_indices.push_back(row_dof);
    offdiag_abs_sum_values.push_back(offdiag_abs_sum);
  }

  predictor.offdiag_abs_sum_raw->insert(offdiag_abs_sum_values, offdiag_abs_sum_indices);
  predictor.offdiag_abs_sum_raw->close();

  mooseAssert(relaxation_parameter > 0.0,
              "The momentum equation relaxation parameter must be positive.");

  for (auto * vector : {predictor.rhs_raw.get(),
                        predictor.pressure_coupling_diagonal_raw.get(),
                        predictor.h_diagonal_correction_raw.get(),
                        predictor.interior_diagonal_raw.get(),
                        predictor.boundary_component_diagonal_raw.get(),
                        predictor.boundary_relax_diagonal_raw.get(),
                        predictor.boundary_a_diagonal_raw.get(),
                        predictor.boundary_dominance_diagonal_raw.get(),
                        predictor.offdiag_abs_sum_raw.get()})
    vector->close();
  predictor.relaxation_parameter = relaxation_parameter;

  // The components are assembled in order, so the last one completes the vector operator.
  if (system_i + 1 == _momentum_systems.size())
    finalizeCachedMomentumPredictorOperators();
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::finalizeCachedMomentumPredictorOperators()
{
  mooseAssert(_momentum_predictor_operators.size() == _momentum_systems.size() &&
                  std::all_of(_momentum_predictor_operators.begin(),
                              _momentum_predictor_operators.end(),
                              [](const auto & op) { return static_cast<bool>(op); }),
              "Every momentum component must be assembled before the predictor is finalized.");

  struct BoundaryFaceCoefficients
  {
    std::vector<Real> coefficient;
    std::vector<dof_id_type> dof;
    std::vector<bool> has_component;
  };

  const auto n_components = _momentum_predictor_operators.size();
  std::unordered_map<dof_id_type, BoundaryFaceCoefficients> boundary_faces;
  for (const auto system_i : index_range(_momentum_predictor_operators))
  {
    auto & predictor = *_momentum_predictor_operators[system_i];
    for (auto * vector : {predictor.boundary_component_diagonal_raw.get(),
                          predictor.boundary_relax_diagonal_raw.get(),
                          predictor.boundary_a_diagonal_raw.get(),
                          predictor.boundary_dominance_diagonal_raw.get()})
    {
      vector->zero();
      vector->close();
    }

    for (const auto & contribution : predictor.boundary_matrix_contributions)
    {
      predictor.boundary_component_diagonal_raw->add(contribution.dof, contribution.contribution);

      auto & face = boundary_faces[contribution.face_id];
      if (face.coefficient.empty())
      {
        face.coefficient.assign(n_components, 0.0);
        face.dof.assign(n_components, DofObject::invalid_id);
        face.has_component.assign(n_components, false);
      }

      face.coefficient[system_i] += contribution.contribution;
      face.dof[system_i] = contribution.dof;
      face.has_component[system_i] = true;
    }
  }

  for (const auto & [face_id, face] : boundary_faces)
  {
    (void)face_id;
    Real component_sum = 0.0;
    Real max_component_magnitude = 0.0;
    Real min_component = std::numeric_limits<Real>::max();
    unsigned int n_present_components = 0;

    for (const auto component_i : index_range(face.coefficient))
      if (face.has_component[component_i])
      {
        const Real coefficient = face.coefficient[component_i];
        component_sum += coefficient;
        max_component_magnitude = std::max(max_component_magnitude, std::abs(coefficient));
        min_component = std::min(min_component, coefficient);
        ++n_present_components;
      }

    if (!n_present_components)
      continue;

    const Real component_average = component_sum / n_present_components;
    for (const auto component_i : index_range(face.coefficient))
      if (face.has_component[component_i])
      {
        const auto dof = face.dof[component_i];
        auto & predictor = *_momentum_predictor_operators[component_i];
        predictor.boundary_dominance_diagonal_raw->add(dof, max_component_magnitude);
        predictor.boundary_relax_diagonal_raw->add(dof, min_component);
        predictor.boundary_a_diagonal_raw->add(dof, component_average);
      }
  }

  for (const auto system_i : index_range(_momentum_predictor_operators))
  {
    auto & predictor = *_momentum_predictor_operators[system_i];

    for (auto * vector : {predictor.boundary_component_diagonal_raw.get(),
                          predictor.boundary_relax_diagonal_raw.get(),
                          predictor.boundary_a_diagonal_raw.get(),
                          predictor.boundary_dominance_diagonal_raw.get()})
      vector->close();

    const auto first_local = predictor.interior_diagonal_raw->first_local_index();
    const auto last_local = predictor.interior_diagonal_raw->last_local_index();
    std::vector<dof_id_type> indices;
    std::vector<Real> pressure_diagonal;
    std::vector<Real> h_diagonal_correction;
    indices.reserve(last_local - first_local);
    pressure_diagonal.reserve(last_local - first_local);
    h_diagonal_correction.reserve(last_local - first_local);

    {
      PetscVectorReader interior_diagonal_reader(*predictor.interior_diagonal_raw);
      PetscVectorReader boundary_component_diagonal_reader(
          *predictor.boundary_component_diagonal_raw);
      PetscVectorReader boundary_relax_diagonal_reader(*predictor.boundary_relax_diagonal_raw);
      PetscVectorReader boundary_a_diagonal_reader(*predictor.boundary_a_diagonal_raw);
      PetscVectorReader boundary_dominance_diagonal_reader(
          *predictor.boundary_dominance_diagonal_raw);
      PetscVectorReader offdiag_abs_sum_reader(*predictor.offdiag_abs_sum_raw);

      for (const auto dof : make_range(first_local, last_local))
      {
        const Real interior_diagonal = interior_diagonal_reader(dof);
        const Real boundary_a_diagonal = boundary_a_diagonal_reader(dof);
        const Real diagonal_for_dominance =
            interior_diagonal + boundary_dominance_diagonal_reader(dof);

        Real predictor_diagonal =
            std::max(std::abs(diagonal_for_dominance), offdiag_abs_sum_reader(dof));

        predictor_diagonal /= predictor.relaxation_parameter;
        predictor_diagonal -= boundary_relax_diagonal_reader(dof);

        indices.push_back(dof);
        pressure_diagonal.push_back(predictor_diagonal + boundary_a_diagonal);
        h_diagonal_correction.push_back(interior_diagonal - predictor_diagonal -
                                        boundary_component_diagonal_reader(dof) +
                                        boundary_a_diagonal);
      }
    }

    predictor.pressure_coupling_diagonal_raw->zero();
    predictor.pressure_coupling_diagonal_raw->close();
    predictor.pressure_coupling_diagonal_raw->insert(pressure_diagonal, indices);
    predictor.pressure_coupling_diagonal_raw->close();

    predictor.h_diagonal_correction_raw->zero();
    predictor.h_diagonal_correction_raw->close();
    predictor.h_diagonal_correction_raw->insert(h_diagonal_correction, indices);
    predictor.h_diagonal_correction_raw->close();
  }

  _predictor_ready = true;
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::computePredictorOperator(const unsigned int system_i,
                                                                       NumericVector<Number> & HbyA,
                                                                       NumericVector<Number> & Ainv)
{
  mooseAssert(_predictor_ready,
              "The momentum predictor operator must be assembled and finalized before H/A.");
  const auto * const predictor = cachedMomentumPredictorOperator(system_i);
  mooseAssert(predictor, "The momentum predictor operator is required before computing HbyA.");

  computePredictorOperatorBase(system_i, HbyA, Ainv, predictor->rhs_raw.get());

  // fvMatrix::H(): the component boundary diagonal is swapped for the averaged one.
  NumericVector<Number> & current_local_solution =
      *(_momentum_implicit_systems[system_i]->current_local_solution);
  auto h_correction = current_local_solution.zero_clone();
  auto * h_correction_petsc = dynamic_cast<PetscVector<Number> *>(h_correction.get());
  mooseAssert(h_correction_petsc,
              "The vectors used in the RhieChowMassFlux objects need to be convertable "
              "to PetscVectors!");
  h_correction_petsc->pointwise_mult(*predictor->h_diagonal_correction_raw, current_local_solution);
  h_correction_petsc->close();
  HbyA.add(1.0, *h_correction);

  // fvMatrix::A(): the relaxed interior diagonal with the averaged boundary diagonal.
  Ainv = *predictor->pressure_coupling_diagonal_raw;
  Ainv.close();
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::maxCourant(const Real dt) const
{
  if (dt <= 0.0)
    return 0.0;

  std::unordered_map<dof_id_type, Real> cell_flux_sum;
  cell_flux_sum.reserve(_flow_face_info.size());

  for (const auto * fi : _flow_face_info)
  {
    if (!fi)
      continue;

    const Real face_measure = NS::faceMeasure(*fi);
    const Real face_volumetric_flux = std::abs(getVolumetricFaceFlux(*fi));
    const Real volumetric_flux = face_volumetric_flux * face_measure;

    if (fi->elemPtr() && hasBlocks(fi->elemPtr()->subdomain_id()))
      cell_flux_sum[fi->elemPtr()->id()] += volumetric_flux;

    if (fi->neighborPtr() && hasBlocks(fi->neighborPtr()->subdomain_id()))
      cell_flux_sum[fi->neighborPtr()->id()] += volumetric_flux;
  }

  Real max_courant = 0.0;
  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
  {
    if (!hasBlocks(elem_info->subdomain_id()))
      continue;

    const Real volume = elem_info->volume() * elem_info->coordFactor();
    if (volume <= libMesh::TOLERANCE)
      continue;

    const auto it = cell_flux_sum.find(elem_info->elem()->id());
    if (it == cell_flux_sum.end())
      continue;

    const Real courant = 0.5 * dt * it->second / volume;
    max_courant = std::max(max_courant, courant);
  }

  // The subcycle count it sizes must agree on every rank
  comm().max(max_courant);
  return max_courant;
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::updatePressureBoundaryNormalGradients()
{
  mooseAssert(_pressure_predictor_face_state_valid,
              "The predictor face state must be published before the boundary gradients.");

  const auto time_arg = Moose::currentState();
  Real max_boundary_normal_ainv = 0.0;

  for (const auto * fi : _flow_face_info)
    if (!_vel[0]->isInternalFace(*fi))
      max_boundary_normal_ainv =
          std::max(max_boundary_normal_ainv, std::abs(pressureBoundaryNormalAinv(fi)));

  // Treat a boundary pressure coefficient as degenerate only relative to the scale of the
  // current boundary operator. An absolute libMesh::TOLERANCE cutoff is too aggressive for
  // high-density-ratio outlet faces, where valid normal Ainv values can be smaller than 1e-6.
  const Real degenerate_normal_ainv_tol =
      std::max(std::numeric_limits<Real>::min(),
               std::numeric_limits<Real>::epsilon() * std::max(1.0, max_boundary_normal_ainv));

  for (const auto * fi : _flow_face_info)
  {
    if (_vel[0]->isInternalFace(*fi))
      continue;

    const Real normal_ainv = pressureBoundaryNormalAinv(fi);
    if (!std::isfinite(normal_ainv) || std::abs(normal_ainv) <= degenerate_normal_ainv_tol)
    {
      _pressure_boundary_normal_gradient[fi->id()] = 0.0;
      continue;
    }

    const Real required_pressure_flux =
        _phiHbyA_flux.evaluate(fi) + pressureBoundaryTargetFlux(fi, time_arg);
    // The diffuse pressure equation stores area-integrated phi; the gradient wants a density.
    const Real face_measure = NS::faceMeasure(*fi);
    if (face_measure <= libMesh::TOLERANCE)
    {
      _pressure_boundary_normal_gradient[fi->id()] = 0.0;
      continue;
    }

    const Real required_pressure_flux_density = required_pressure_flux / face_measure;
    _pressure_boundary_normal_gradient[fi->id()] = -required_pressure_flux_density / normal_ainv;
  }

  _pressure_boundary_normal_gradient_valid = true;
  ++_face_mass_flux_generation;
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::updateVelocityBoundaryState()
{
  const auto time_arg = Moose::currentState();

  _velocity_boundary_state_valid = false;
  for (auto & component_face_values : _boundary_velocity_face_values)
    component_face_values.clear();

  for (const auto * fi : _flow_face_info)
  {
    if (_vel[0]->isInternalFace(*fi))
      continue;

    RealVectorValue density_times_velocity;
    const Real boundary_normal_multiplier = boundaryNormalMultiplier(fi);
    const Moose::FaceArg boundary_face = makeCenteredFaceArg(fi, boundaryElemInfo(fi).elem());
    const Real face_rho = _rho(boundary_face, time_arg);

    for (const auto component : index_range(_vel))
    {
      if (auto * bc_pointer = velocityBoundaryCondition(fi, component))
      {
        _boundary_velocity_face_values[component][fi->id()] = bc_pointer->computeBoundaryValue();
      }
      else
      {
        const ElemInfo & elem_info = boundaryElemInfo(fi);
        _boundary_velocity_face_values[component][fi->id()] =
            _vel[component]->getElemValue(elem_info, time_arg);
      }

      density_times_velocity(component) = boundary_normal_multiplier * face_rho *
                                          _boundary_velocity_face_values[component][fi->id()];
    }

    // Before a pressure correction has been solved, publish the predictor-side
    // boundary phi from the refreshed patch state. But once the pressure
    // equation has produced a corrected face flux, preserve that phi so the
    // boundary face correction stays matched to the cell writeback instead of
    // being overwritten by a later cache refresh.
    if (!_pressure_equation_flux_valid)
      _face_mass_flux[fi->id()] = density_times_velocity * fi->normal();
  }

  _velocity_boundary_state_valid = true;
  cacheCurrentCorrectedVolumetricFlux();
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::boundaryVelocityValue(
    const FaceInfo * fi, const unsigned int component, const Moose::StateArg & time_arg) const
{
  if (_velocity_boundary_state_valid)
    if (const auto it = _boundary_velocity_face_values[component].find(fi->id());
        it != _boundary_velocity_face_values[component].end())
      return it->second;

  return RhieChowMassFlux::boundaryVelocityValue(fi, component, time_arg);
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::boundaryVolumetricFluxTarget(
    const FaceInfo * fi, const Moose::StateArg & time_arg) const
{
  mooseAssert(fi && !_vel[0]->isInternalFace(*fi),
              "boundaryVolumetricFluxTarget should only be called on boundary faces.");

  const Real boundary_normal_multiplier = boundaryNormalMultiplier(fi);

  RealVectorValue face_velocity;
  for (const auto component : index_range(_vel))
    face_velocity(component) =
        boundary_normal_multiplier * boundaryVelocityValue(fi, component, time_arg);

  return face_velocity * fi->normal();
}

RealVectorValue
ConservativeDiffuseInterfaceRhieChowMassFlux::interpolateMomentumVectorReadersToFace(
    const FaceInfo * fi, const std::vector<PetscVectorReader> & readers) const
{
  using namespace Moose::FV;

  RealVectorValue face_value;

  if (!fi || readers.size() < _dim)
    return face_value;

  if (_vel[0]->isInternalFace(*fi))
  {
    const auto & elem_info = *fi->elemInfo();
    const auto & neighbor_info = *fi->neighborInfo();

    for (const auto dim_i : make_range(_dim))
    {
      const auto elem_dof = elem_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
      const auto neighbor_dof =
          neighbor_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
      face_value(dim_i) =
          linearInterpolation(readers[dim_i](elem_dof), readers[dim_i](neighbor_dof), *fi, true);
    }

    return face_value;
  }

  const ElemInfo & elem_info = boundaryElemInfo(fi);
  for (const auto dim_i : make_range(_dim))
  {
    const auto elem_dof = elem_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
    face_value(dim_i) = readers[dim_i](elem_dof);
  }

  return face_value;
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::computeFaceNormalAinv(
    const RealVectorValue & face_ainv, const RealVectorValue & face_normal) const
{
  Real normal_ainv = 0.0;
  for (const auto dim_i : make_range(_dim))
    normal_ainv += face_ainv(dim_i) * face_normal(dim_i) * face_normal(dim_i);

  return normal_ainv;
}

InputParameters
ConservativeDiffuseInterfaceRhieChowMassFlux::validParams()
{
  InputParameters params = RhieChowMassFlux::validParams();

  params.addClassDescription(
      "Rhie-Chow face-flux provider for the reduced-pressure diffuse-interface two-phase solve: "
      "the momentum predictor is assembled without the pressure, gravity and capillary forces, "
      "which enter the pressure equation as face fluxes, and the face mass flux is the one the "
      "limited volume-fraction transport advanced with.");

  params.addParam<MooseFunctorName>(
      "vof_rho_phi_functor",
      "rho_phi",
      "Optional alpha-consistent density flux published by the diffuse-interface VOF transport. "
      "When present, downstream advection queries use this face mass flux instead of the raw "
      "Rhie-Chow density interpolation.");
  params.addParam<RealVectorValue>(
      "gravity",
      RealVectorValue(),
      "Gravity vector used by the reduced-pressure hydrostatic pressure initializer.");
  params.addParam<Point>(
      "reference_pressure_point",
      Point(0, 0, 0),
      "Reference point used to compute the reduced-pressure head gh = g.(x-x_ref).");
  params.addParam<MooseFunctorName>(
      "surface_tension_force_density_functor",
      "",
      "Optional surface tension force density functor sigma*kappa*grad(alpha), as published by "
      "ConservativeDiffuseInterfaceSurfaceTensionFunctorMaterial. When supplied, the "
      "diffuse-interface pressure coupling carries the capillary force.");
  params.set<MooseEnum>("pressure_projection_method") =
      MooseEnum("standard consistent", "consistent");

  return params;
}

ConservativeDiffuseInterfaceRhieChowMassFlux::ConservativeDiffuseInterfaceRhieChowMassFlux(
    const InputParameters & params)
  : RhieChowMassFlux(params),
    _phiHbyA_flux(_moose_mesh,
                  blockIDs(),
                  "phiHbyA_flux",
                  [this](const FaceInfo & fi)
                  {
                    // Total by construction: a face the predictor has not reached reads as zero,
                    // which is what the stored field gave once zero-initialised.
                    const auto at = [&fi](const auto & field)
                    {
                      const auto it = field.find(fi.id());
                      return it == field.end() ? 0.0 : it->second;
                    };
                    return at(_pressure_predictor_base_flux) - at(_phig_flux);
                  }),
    _phig_flux(_moose_mesh, blockIDs(), "phig_flux"),
    _pressure_Ainv(_moose_mesh, blockIDs(), "pressure_Ainv"),
    _pressure_boundary_normal_gradient(
        _moose_mesh, blockIDs(), "pressure_boundary_normal_gradient"),
    _boundary_velocity_face_values(_dim),
    _corrected_face_phi(_moose_mesh,
                        blockIDs(),
                        "corrected_face_phi",
                        [this](const FaceInfo & fi) { return correctedTransportPhi(fi); }),
    _previous_timestep_corrected_face_phi(
        _moose_mesh, blockIDs(), "previous_timestep_corrected_face_phi"),
    _lagged_corrected_face_phi(_moose_mesh, blockIDs(), "lagged_corrected_face_phi"),
    _pressure_coupled_cell_reconstruction_scalar(
        _moose_mesh, blockIDs(), "pressure_coupled_cell_reconstruction_scalar"),
    _gravity(getParam<RealVectorValue>("gravity")),
    _reference_pressure_point(getParam<Point>("reference_pressure_point")),
    _vof_rho_phi_name(getParam<MooseFunctorName>("vof_rho_phi_functor")),
    _surface_tension_force_name(
        getParam<MooseFunctorName>("surface_tension_force_density_functor")),
    _vof_rho_phi(nullptr)
{
  for (const auto tid : make_range(libMesh::n_threads()))
  {
    UserObject::_subproblem.addFunctor("pressure_Ainv", _pressure_Ainv, tid);
    UserObject::_subproblem.addFunctor("phiHbyA", _phiHbyA_flux, tid);
    UserObject::_subproblem.addFunctor(
        "pressure_boundary_normal_gradient", _pressure_boundary_normal_gradient, tid);
    UserObject::_subproblem.addFunctor("corrected_face_phi", _corrected_face_phi, tid);
    UserObject::_subproblem.addFunctor(
        "previous_timestep_corrected_face_phi", _previous_timestep_corrected_face_phi, tid);
    UserObject::_subproblem.addFunctor(
        "lagged_corrected_face_phi", _lagged_corrected_face_phi, tid);
    UserObject::_subproblem.addFunctor("pressure_coupled_cell_reconstruction_scalar",
                                       _pressure_coupled_cell_reconstruction_scalar,
                                       tid);
  }

  initializeAdditionalPressureFluxStorage();
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::getMassFlux(const FaceInfo & fi) const
{
  // Deliberately not the inherited _face_mass_flux, which this class also maintains as the
  // restartable rho_f * phi_f. Momentum advection needs the density flux that the limited volume
  // fraction update actually transported, and across a large density jump the two differ: the
  // transported product does not factor into one face density times one face volumetric flux.
  // _face_mass_flux remains what the base class and restart expect; this is what advection needs.
  const Real face_measure = NS::faceMeasure(fi);
  return face_measure > 0.0 ? vofRhoPhiIntegrated(fi) / face_measure : 0.0;
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::commitAcceptedTimestepTransportHistory()
{
  for (const auto * fi : _flow_face_info)
    _previous_timestep_corrected_face_phi[fi->id()] = _corrected_face_phi.evaluate(fi);
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::vofRhoPhiIntegrated(const FaceInfo & fi) const
{
  return evaluateFaceScalarFunctor(_vof_rho_phi, &fi, Moose::currentState());
}

RealVectorValue
ConservativeDiffuseInterfaceRhieChowMassFlux::reducedPressureMomentumPredictorForceDensity(
    const ElemInfo & elem_info, const Moose::StateArg & time_arg) const
{
  if (!hasBlocks(elem_info.subdomain_id()))
    return RealVectorValue();

  RealVectorValue force_density = reconstructCellVectorFromFaceNormalScalars(
      _moose_mesh,
      _dim,
      &elem_info,
      [this, &time_arg](const FaceInfo & fi_loc, const Real, Real & normal_force_density)
      {
        if (!_vel[0]->isInternalFace(fi_loc))
          return false;

        normal_force_density = -computeFaceNormalPressureGradient(&fi_loc, time_arg);
        if (useExplicitHydrostaticPredictorForce())
          normal_force_density +=
              computeDiscreteHydrostaticPredictorFaceNormalForceDensity(&fi_loc, time_arg);

        return true;
      });

  force_density += surfaceTensionForceDensity(elem_info, time_arg);

  return force_density;
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::predictorVelocityComponent(
    const ElemInfo & elem_info, const unsigned int component) const
{
  mooseAssert(component < _dim, "Momentum component index out of range.");

  const auto dof = elem_info.dofIndices()[_global_momentum_system_numbers[component]][0];
  return -(*_HbyA_raw[component])(dof);
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::getVolumetricFaceFlux(const FaceInfo & fi) const
{
  return _corrected_face_phi.evaluate(&fi);
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::initializeAdditionalPressureFluxStorage()
{
  for (const auto * fi : _fe_problem.mesh().faceInfo())
  {
    _previous_timestep_corrected_face_phi[fi->id()] = 0.0;
    _lagged_corrected_face_phi[fi->id()] = 0.0;
    _pressure_coupled_cell_reconstruction_scalar[fi->id()] = 0.0;
    _phig_flux[fi->id()] = 0.0;
    _pressure_predictor_base_flux[fi->id()] = 0.0;
    _pressure_Ainv[fi->id()] = RealVectorValue();
    _pressure_boundary_normal_gradient[fi->id()] = 0.0;
  }

  for (auto & component_face_values : _boundary_velocity_face_values)
    component_face_values.clear();
  _velocity_boundary_state_valid = false;
  _pressure_boundary_normal_gradient_valid = false;
  _pressure_predictor_face_state_valid = false;
  _pressure_coupled_velocity_correction_valid = false;
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::writePressureCorrectedVelocityToMomentumSolution(
    const Moose::StateArg & time_arg)
{
  std::vector<std::unique_ptr<NumericVector<Number>>> pressure_corrected_solution;
  pressure_corrected_solution.reserve(_momentum_implicit_systems.size());

  for (const auto system_i : index_range(_momentum_implicit_systems))
  {
    auto * momentum_system = _momentum_implicit_systems[system_i];
    mooseAssert(momentum_system && momentum_system->current_local_solution,
                "The requested momentum component is not linked to "
                "ConservativeDiffuseInterfaceRhieChowMassFlux.");
    pressure_corrected_solution.push_back(momentum_system->current_local_solution->zero_clone());
  }

  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
  {
    if (!hasBlocks(elem_info->subdomain_id()))
      continue;

    const RealVectorValue pressure_delta =
        reconstructPressureCoupledCellVelocityDelta(elem_info, time_arg);

    for (const auto system_i : index_range(_momentum_implicit_systems))
    {
      const auto & dof_indices = elem_info->dofIndices()[_global_momentum_system_numbers[system_i]];
      if (dof_indices.empty())
        continue;

      pressure_corrected_solution[system_i]->set(dof_indices[0],
                                                 predictorVelocityComponent(*elem_info, system_i) +
                                                     pressure_delta(system_i));
    }
  }

  for (const auto system_i : index_range(_momentum_implicit_systems))
  {
    pressure_corrected_solution[system_i]->close();
    *(_momentum_implicit_systems[system_i]->solution) = *pressure_corrected_solution[system_i];
    _momentum_implicit_systems[system_i]->update();
    _momentum_systems[system_i]->setSolution(
        *_momentum_implicit_systems[system_i]->current_local_solution);
  }
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::transportMassFluxDensityFromVolumetricPhi(
    const FaceInfo * fi, const Real volumetric_phi, const Moose::StateArg & time_arg) const
{
  const Real face_rho = interpolateFaceDensity(fi, time_arg);
  return face_rho * volumetric_phi;
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::pressureBoundaryTargetFlux(
    const FaceInfo * fi, const Moose::StateArg & time_arg) const
{
  if (!fi)
    return 0.0;

  return boundaryVolumetricFluxTarget(fi, time_arg) * NS::faceMeasure(*fi);
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::pressureBoundaryNormalAinv(const FaceInfo * fi) const
{
  if (!fi)
    return 0.0;

  return computeFaceNormalAinv(libmesh_map_find(_pressure_Ainv, fi->id()), fi->normal());
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::pressureFaceNormalAinv(
    const FaceInfo * fi, const Moose::StateArg & time_arg) const
{
  if (!fi)
    return 0.0;

  if (!_vel[0]->isInternalFace(*fi))
    return pressureBoundaryNormalAinv(fi);

  // The pressure writeback must divide by the same face coefficient consumed
  // by the pressure diffusion kernel. Read that coefficient back through the
  // pressure-space face functor rather than projecting the cached raw Ainv map
  // directly, so the writeback stays aligned with the operator contract if the
  // pressure face interpolation changes.
  const RealVectorValue pressure_face_rau = _pressure_Ainv(makeCenteredFaceArg(fi), time_arg);
  return computeFaceNormalAinv(pressure_face_rau, fi->normal());
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::correctedTransportPhi(const FaceInfo & fi) const
{
  // The corrected transport flux is the predictor flux plus whatever the pressure equation added,
  // both of which the provider already holds. Expressing it rather than storing it means a
  // consumer cannot read a value that was correct for an earlier predictor state: that staleness
  // is what let a transient outlet backflow decision persist.
  // A face the flux fields have not reached yet reads as no flux, which is what the stored field
  // used to give once it had been zero-initialised. Keeping the definition total means the early
  // calls during initialisation behave as they did.
  const auto face_flux_of = [&fi](const auto & field)
  {
    const auto it = field.find(fi.id());
    return it == field.end() ? 0.0 : it->second;
  };

  if (!_pressure_predictor_face_state_valid)
  {
    // Before a predictor exists the only mass flux is the one initFaceMassFlux built as
    // face_rho * boundary velocity, formed with a single face density, so dividing it back out
    // recovers that velocity rather than inverting a transported product. This is the one place
    // the division is exact, which is why it is confined to this branch.
    const Real face_rho = predictorFaceDensity(&fi, Moose::currentState());
    return std::abs(face_rho) > libMesh::TOLERANCE ? face_flux_of(_face_mass_flux) / face_rho : 0.0;
  }

  const Real face_measure = NS::faceMeasure(fi);
  if (face_measure <= libMesh::TOLERANCE)
    return 0.0;

  // The diffuse pressure equation stores area-integrated phi; transport wants a flux density.
  const Real predictor_phi = -_phiHbyA_flux.evaluate(&fi) / face_measure;
  const Real pressure_equation_phi =
      _pressure_equation_flux_valid ? face_flux_of(_pressure_equation_flux) / face_measure : 0.0;
  return predictor_phi + pressure_equation_phi;
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::writebackReconstructionScalar(
    const FaceInfo * fi,
    const Moose::StateArg & time_arg,
    const Real degenerate_normal_pressure_ainv_tol) const
{
  if (!fi || !_pressure_predictor_face_state_valid || !_pressure_equation_flux_valid)
    return 0.0;

  // The diffuse pressure equation stores area-integrated phi. Convert back to the
  // normal-flux-density transport convention only at the writeback boundary.
  const Real face_measure = NS::faceMeasure(*fi);
  if (face_measure <= libMesh::TOLERANCE)
    return 0.0;

  const Real pressure_writeback_phi =
      libmesh_map_find(_pressure_equation_flux, fi->id()) / face_measure +
      libmesh_map_find(_phig_flux, fi->id()) / face_measure;
  const Real normal_pressure_ainv = pressureFaceNormalAinv(fi, time_arg);

  return std::isfinite(normal_pressure_ainv) &&
                 std::abs(normal_pressure_ainv) > degenerate_normal_pressure_ainv_tol
             ? pressure_writeback_phi / normal_pressure_ainv
             : 0.0;
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::maxPressureFaceNormalAinv(
    const Moose::StateArg & time_arg) const
{
  Real max_abs_normal_pressure_ainv = 0.0;
  for (const auto * fi : _flow_face_info)
  {
    const Real normal_pressure_ainv = pressureFaceNormalAinv(fi, time_arg);
    if (std::isfinite(normal_pressure_ainv))
      max_abs_normal_pressure_ainv =
          std::max(max_abs_normal_pressure_ainv, std::abs(normal_pressure_ainv));
  }

  return max_abs_normal_pressure_ainv;
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::cacheCurrentCorrectedVolumetricFlux(
    const Real degenerate_normal_pressure_ainv_tol)
{
  const auto time_arg = Moose::currentState();
  // The corrected face flux stays authoritative for continuity and transport; this scalar only
  // feeds the cell-centered velocity writeback and is never projected back onto the faces.
  for (const auto * fi : _flow_face_info)
    _pressure_coupled_cell_reconstruction_scalar[fi->id()] =
        writebackReconstructionScalar(fi, time_arg, degenerate_normal_pressure_ainv_tol);
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::freezeLaggedCorrectedFaceFlux(
    const bool use_previous_timestep_flux)
{
  for (const auto * fi : _flow_face_info)
    _lagged_corrected_face_phi[fi->id()] =
        use_previous_timestep_flux
            ? libmesh_map_find(_previous_timestep_corrected_face_phi, fi->id())
            : _corrected_face_phi.evaluate(fi);
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::meshChanged()
{
  RhieChowMassFlux::meshChanged();
  _momentum_predictor_operators.clear();
  _predictor_ready = false;
  initializeAdditionalPressureFluxStorage();
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::initialSetup()
{
  const_cast<MooseLinearVariableFVReal *>(_p)->computeCellGradients();
  RhieChowMassFlux::initialSetup();
  initializeAdditionalPressureFluxStorage();
  if (!_vof_rho_phi_name.empty() &&
      UserObject::_subproblem.hasFunctorWithType<Real>(_vof_rho_phi_name, _tid))
    _vof_rho_phi = &getFunctor<Real>(_vof_rho_phi_name);

  if (!_vof_rho_phi)
    paramError("vof_rho_phi_functor",
               "The conservative diffuse-interface path requires the VOF-owned rhoPhi functor.");

  if (!_surface_tension_force_name.empty())
    _surface_tension_force = &getFunctor<RealVectorValue>(_surface_tension_force_name);
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::initialize()
{
  RhieChowMassFlux::initialize();

  initializeAdditionalPressureFluxStorage();
  commitAcceptedTimestepTransportHistory();
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::cachePressureEquationFlux()
{
  const_cast<MooseLinearVariableFVReal *>(_p)->computeCellGradients();

  for (auto & fi : _flow_face_info)
    _pressure_equation_flux[fi->id()] = computeDiscretePressureFaceFlux(fi);

  _pressure_equation_flux_valid = true;

  const auto time_arg = Moose::currentState();
  const Real max_abs_normal_pressure_ainv = maxPressureFaceNormalAinv(time_arg);
  const Real degenerate_normal_pressure_ainv_tol =
      std::max(std::numeric_limits<Real>::min(),
               pressure_writeback_face_ainv_relative_tolerance * max_abs_normal_pressure_ainv);
  cacheCurrentCorrectedVolumetricFlux(degenerate_normal_pressure_ainv_tol);
  _pressure_coupled_velocity_correction_valid = false;
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::computeHbyA(const bool verbose)
{
  // A new predictor invalidates everything derived from the previous one.
  _pressure_predictor_face_state_valid = false;
  _pressure_boundary_normal_gradient_valid = false;
  RhieChowMassFlux::computeHbyA(verbose);
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::computeFaceMassFlux()
{
  mooseAssert(_pressure_predictor_face_state_valid,
              "The predictor face state must be published before the face mass flux.");

  if (!_pressure_equation_flux_valid)
    cachePressureEquationFlux();

  cacheCurrentCorrectedVolumetricFlux();

  const auto time_arg = Moose::currentState();
  for (const auto * fi : _flow_face_info)
  {
    // Keep the pressure-corrected face flux as the source of truth. Do not reconstruct cell U and
    // project it back to faces here; that round trip is not an exact inverse.
    _face_mass_flux[fi->id()] =
        transportMassFluxDensityFromVolumetricPhi(fi, _corrected_face_phi.evaluate(fi), time_arg);
  }
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::interpolateFaceDensity(
    const FaceInfo * fi, const Moose::StateArg & time_arg) const
{
  if (_vel[0]->isInternalFace(*fi))
    return evaluateFaceScalarFunctor(&_rho, fi, time_arg);

  const ElemInfo & elem_info = boundaryElemInfo(fi);

  if (useConstrainedBoundaryPredictorState(fi))
    return _rho(makeCenteredFaceArg(fi, elem_info.elem()), time_arg);

  return _rho(makeElemArg(elem_info.elem()), time_arg);
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::predictorFaceDensity(
    const FaceInfo * fi, const Moose::StateArg & time_arg) const
{
  // In the diffuse-interface VOF path, rhoPhi is assembled from the primary volumetric
  // flux plus the limited alpha transport flux. It is therefore not generally equal
  // to rho_f * phi, and recovering a predictor "face density" from rhoPhi / phi
  // injects a non-physical variable-density mismatch into the pressure predictor and
  // velocity writeback chain.
  return interpolateFaceDensity(fi, time_arg);
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::facePhysicalVelocityComponent(
    const FaceInfo * fi, const unsigned int component, const Moose::StateArg & time_arg) const
{
  mooseAssert(fi, "FaceInfo should not be null when evaluating a face velocity component.");
  mooseAssert(component < _vel.size(), "Velocity component index out of range.");

  if (_vel[component]->isInternalFace(*fi))
  {
    const auto & elem_info = *fi->elemInfo();
    const auto & neighbor_info = *fi->neighborInfo();
    return Moose::FV::linearInterpolation(_vel[component]->getElemValue(elem_info, time_arg),
                                          _vel[component]->getElemValue(neighbor_info, time_arg),
                                          *fi,
                                          true);
  }

  return boundaryVelocityValue(fi, component, time_arg);
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::interpolatedPhysicalFaceFlux(
    const FaceInfo * fi, const Moose::StateArg & time_arg) const
{
  if (!fi)
    return 0.0;

  RealVectorValue face_velocity;
  for (const auto dim_i : make_range(_dim))
    face_velocity(dim_i) = facePhysicalVelocityComponent(fi, dim_i, time_arg);

  return face_velocity * fi->normal();
}

ConservativeDiffuseInterfaceRhieChowMassFlux::DiffuseFaceOperatorState
ConservativeDiffuseInterfaceRhieChowMassFlux::buildDiffuseFaceOperatorState(
    const FaceInfo * fi,
    const Moose::StateArg & time_arg,
    const std::vector<PetscVectorReader> & raw_ainv_readers) const
{
  DiffuseFaceOperatorState state;
  state.face_normal = fi->normal();
  state.face_rho = interpolateFaceDensity(fi, time_arg);
  state.face_raw_ainv = interpolateMomentumVectorReadersToFace(fi, raw_ainv_readers);
  state.normal_raw_ainv = computeFaceNormalAinv(state.face_raw_ainv, state.face_normal);
  return state;
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::evaluateFaceScalarFunctor(
    const Moose::Functor<Real> * functor,
    const FaceInfo * fi,
    const Moose::StateArg & time_arg) const
{
  if (!functor)
    return 0.0;

  return MetaPhysicL::raw_value((*functor)(makeCenteredFaceArg(fi), time_arg));
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::computeDefaultTransientProjectionVolumetricFlux(
    const FaceInfo * fi, const Moose::StateArg & time_arg, const DiffuseFaceOperatorState & state)
{
  if (!fi || time_arg.state != Moose::currentState().state)
    return 0.0;

  const Real dt = _fe_problem.dt();
  if (dt <= libMesh::TOLERANCE)
    return 0.0;

  if (!_vel[0]->isInternalFace(*fi))
    return 0.0;

  const Real old_corrected_phi = libmesh_map_find(_previous_timestep_corrected_face_phi, fi->id());
  const Real old_interpolated_face_phi = interpolatedPhysicalFaceFlux(fi, Moose::oldState());
  const Real phi_corr = old_corrected_phi - old_interpolated_face_phi;
  const Real phi_scale = std::abs(old_corrected_phi) + std::numeric_limits<Real>::epsilon();
  const Real ddt_coupling_coeff =
      std::max(0.0, 1.0 - std::min(std::abs(phi_corr) / phi_scale, 1.0));

  // Euler ddtCorr contribution: fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), phiCorr) * phiCorr / dt.
  // The pressure corrector adds interpolate(rho*rAU()) times that value to phiHbyA.
  return ddt_coupling_coeff * state.face_rho * state.normal_raw_ainv * phi_corr / dt;
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::computeFaceNormalDensityGradient(
    const FaceInfo * fi, const Moose::StateArg & time_arg) const
{
  if (!fi)
    return 0.0;

  const bool internal = _vel[0]->isInternalFace(*fi);
  return NS::faceNormalScalarGradient(
      _rho, *fi, internal ? nullptr : boundaryElemInfo(fi).elem(), internal, time_arg);
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::computeFaceNormalPressureGradient(
    const FaceInfo * fi, const Moose::StateArg & time_arg) const
{
  if (!fi)
    return 0.0;

  if (_vel[0]->isInternalFace(*fi))
  {
    const Real elem_p = _p->getElemValue(*fi->elemInfo(), time_arg);
    const Real neighbor_p = _p->getElemValue(*fi->neighborInfo(), time_arg);
    const Real normal_spacing = std::abs(fi->dCN() * fi->normal());

    return normal_spacing > libMesh::TOLERANCE ? (neighbor_p - elem_p) / normal_spacing : 0.0;
  }

  if (_pressure_boundary_normal_gradient_valid)
    return libmesh_map_find(_pressure_boundary_normal_gradient, fi->id());

  const ElemInfo & elem_info = boundaryElemInfo(fi);
  const Point fluid_centroid = elem_info.centroid();
  const Point cell_to_face = fi->faceCentroid() - fluid_centroid;
  const Real normal_spacing = std::abs(cell_to_face * fi->normal());

  if (normal_spacing <= libMesh::TOLERANCE)
    return 0.0;

  const Real elem_p = _p->getElemValue(elem_info, time_arg);

  if (auto * bc_pointer = pressureBoundaryCondition(fi))
  {
    const Real bc_sn_grad = bc_pointer->computeBoundaryNormalGradient();
    if (std::isfinite(bc_sn_grad))
      return bc_sn_grad;

    const Real bc_value = bc_pointer->computeBoundaryValue();
    if (std::isfinite(bc_value))
      return (bc_value - elem_p) / normal_spacing;
  }

  return 0.0;
}

bool
ConservativeDiffuseInterfaceRhieChowMassFlux::useExplicitHydrostaticPredictorForce() const
{
  return !_startup_projection_active;
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::
    computeDiscreteHydrostaticPredictorFaceNormalForceDensity(
        const FaceInfo * fi, const Moose::StateArg & time_arg) const
{
  if (!fi || !useExplicitHydrostaticPredictorForce())
    return 0.0;

  const Real ghf = _gravity * (fi->faceCentroid() - _reference_pressure_point);
  return -ghf * computeFaceNormalDensityGradient(fi, time_arg);
}

bool
ConservativeDiffuseInterfaceRhieChowMassFlux::surfaceTensionEnabled() const
{
  return _surface_tension_force;
}

RealVectorValue
ConservativeDiffuseInterfaceRhieChowMassFlux::surfaceTensionForceDensity(
    const ElemInfo & elem_info, const Moose::StateArg & time_arg) const
{
  // The force itself belongs to a functor material. What belongs here is only whether this
  // corrector may use it, which is startup projection state no material can see.
  if (!surfaceTensionEnabled() || _startup_projection_active || !elem_info.elem())
    return RealVectorValue();

  return (*_surface_tension_force)(makeElemArg(elem_info.elem()), time_arg);
}

Real
ConservativeDiffuseInterfaceRhieChowMassFlux::surfaceTensionFaceNormalForceDensity(
    const FaceInfo * fi, const Moose::StateArg & time_arg) const
{
  if (!fi || !surfaceTensionEnabled() || _startup_projection_active)
    return 0.0;

  // An interior face averages the two cell forces; a boundary face names the side that owns the
  // cell quantity, since there is no second cell to average with.
  const bool internal = _vel[0]->isInternalFace(*fi);
  const Moose::FaceArg face_arg =
      makeCenteredFaceArg(fi, internal ? nullptr : boundaryElemInfo(fi).elem());
  return (*_surface_tension_force)(face_arg, time_arg) * fi->normal();
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::updateAdditionalPressureFluxFunctors(
    const bool /*verbose*/)
{
  _pressure_coupled_velocity_correction_valid = false;
  _pressure_predictor_face_state_valid = false;

  const auto time_arg = Moose::currentState();

  mooseAssert(_Ainv_raw.size() >= _dim, "computeHbyA must run before the face operators.");
  std::vector<PetscVectorReader> raw_ainv_readers;
  raw_ainv_readers.reserve(_dim);
  for (const auto dim_i : make_range(_dim))
    raw_ainv_readers.emplace_back(*_Ainv_raw[dim_i]);

  std::vector<PetscVectorReader> hbya_readers;
  hbya_readers.reserve(_HbyA_raw.size());
  for (const auto & hbya_raw : _HbyA_raw)
    hbya_readers.emplace_back(*hbya_raw);

  for (const auto * fi : _flow_face_info)
  {
    const auto predictor_state =
        buildMomentumPredictorFaceState(fi, hbya_readers, nullptr, time_arg);

    const Real face_measure = NS::faceMeasure(*fi);
    // FaceMomentumPredictorState carries H/A with the sign it is stored in, so that the velocity
    // u = -(H/A) - Ainv grad(p). The reduced-pressure predictor flux is a velocity, and the
    // bookkeeping below negates it once more when forming the base flux, so convert here rather
    // than asking the shared producer to hold a second sign convention for this consumer alone.
    const Real predictor_operator_volumetric_flux = -(predictor_state.face_hbya * fi->normal());
    const auto state = buildDiffuseFaceOperatorState(fi, time_arg, raw_ainv_readers);
    const RealVectorValue pressure_face_raw_ainv = state.face_raw_ainv;
    const Real normal_pressure_ainv =
        computeFaceNormalAinv(pressure_face_raw_ainv, state.face_normal);
    Real transient_projection_flux = 0.0;
    Real hydrostatic_flux = 0.0;
    Real surface_tension_flux = 0.0;

    if (!_startup_projection_active)
    {
      transient_projection_flux =
          -computeDefaultTransientProjectionVolumetricFlux(fi, time_arg, state);

      const Real ghf = _gravity * (fi->faceCentroid() - _reference_pressure_point);
      const Real sn_grad_rho = computeFaceNormalDensityGradient(fi, time_arg);
      hydrostatic_flux = -ghf * sn_grad_rho * normal_pressure_ainv;

      surface_tension_flux =
          surfaceTensionFaceNormalForceDensity(fi, time_arg) * normal_pressure_ainv;
    }

    const Real phig_flux = (hydrostatic_flux + surface_tension_flux) * face_measure;
    const Real pressure_predictor_base_flux =
        (predictor_operator_volumetric_flux + transient_projection_flux) * face_measure;

    _phig_flux[fi->id()] = phig_flux;
    _pressure_Ainv[fi->id()] = pressure_face_raw_ainv;
    _pressure_predictor_base_flux[fi->id()] = -pressure_predictor_base_flux;
  }

  _pressure_predictor_face_state_valid = true;
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::updatePressureCoupledVelocityCorrectionFaceField(
    const Moose::StateArg & time_arg)
{
  if (!_pressure_equation_flux_valid)
    cachePressureEquationFlux();

  // The tolerance is documented as a fraction of the active-face maximum.
  // Using max(1, max_abs_normal_pressure_ainv) turns it into an unintended
  // absolute cutoff whenever the live face-normal Ainv values are << 1, which
  // suppresses writeback reconstruction on the very interface faces we need to
  // correct. Keep the cutoff relative to the active maximum instead.
  const Real max_abs_normal_pressure_ainv = maxPressureFaceNormalAinv(time_arg);
  const Real degenerate_normal_pressure_ainv_tol =
      std::max(std::numeric_limits<Real>::min(),
               pressure_writeback_face_ainv_relative_tolerance * max_abs_normal_pressure_ainv);

  cacheCurrentCorrectedVolumetricFlux(degenerate_normal_pressure_ainv_tol);

  _pressure_coupled_velocity_correction_valid = true;
}

bool
ConservativeDiffuseInterfaceRhieChowMassFlux::useConstrainedBoundaryPredictorState(
    const FaceInfo * fi) const
{
  if (!fi || _vel[0]->isInternalFace(*fi))
    return false;

  bool use_constrained_boundary_state = _vel[0]->isDirichletBoundaryFace(*fi);
  if (!use_constrained_boundary_state)
  {
    for (const auto dim_i : make_range(_dim))
      if (auto * bc_pointer = velocityBoundaryCondition(fi, dim_i))
      {
        if (dynamic_cast<LinearFVPressureInletOutletVelocityBC *>(bc_pointer))
          continue;

        if (auto * adv_diff_bc = dynamic_cast<LinearFVAdvectionDiffusionBC *>(bc_pointer))
        {
          if (adv_diff_bc->computeBoundaryGradientMatrixContribution() > 0.0)
          {
            use_constrained_boundary_state = true;
            break;
          }
        }
      }
  }

  return use_constrained_boundary_state;
}

RealVectorValue
ConservativeDiffuseInterfaceRhieChowMassFlux::reconstructFaceNormalScalarToCellVector(
    const ElemInfo * elem_info,
    const Moose::StateArg & time_arg,
    const FaceScalarField & scalar_field) const
{
  return reconstructCellVectorFromFaceNormalScalars(
      _moose_mesh,
      _dim,
      elem_info,
      [this, &scalar_field, &time_arg](const FaceInfo & fi_loc, const Real, Real & normal_scalar)
      {
        normal_scalar = scalar_field(makeCenteredFaceArg(&fi_loc), time_arg);
        return true;
      });
}

RealVectorValue
ConservativeDiffuseInterfaceRhieChowMassFlux::reconstructBasePressureCoupledCellVelocityDelta(
    const ElemInfo * elem_info, const Moose::StateArg & time_arg) const
{
  if (!elem_info || !_pressure_coupled_velocity_correction_valid)
    return RealVectorValue();

  const RealVectorValue reconstructed_operator = reconstructFaceNormalScalarToCellVector(
      elem_info, time_arg, _pressure_coupled_cell_reconstruction_scalar);

  RealVectorValue pressure_delta;
  for (const auto system_i : index_range(_momentum_implicit_systems))
  {
    const auto & dof_indices = elem_info->dofIndices()[_global_momentum_system_numbers[system_i]];
    if (dof_indices.empty())
      continue;

    pressure_delta(system_i) =
        (*_Ainv_raw[system_i])(dof_indices[0]) * reconstructed_operator(system_i);
  }

  return pressure_delta;
}

RealVectorValue
ConservativeDiffuseInterfaceRhieChowMassFlux::reconstructPressureCoupledCellVelocityDelta(
    const ElemInfo * elem_info, const Moose::StateArg & time_arg) const
{
  if (!elem_info || !_pressure_coupled_velocity_correction_valid)
    return RealVectorValue();

  return reconstructBasePressureCoupledCellVelocityDelta(elem_info, time_arg);
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::computePressureCorrectedCellVelocity()
{
  const auto time_arg = Moose::currentState();
  updatePressureCoupledVelocityCorrectionFaceField(time_arg);
  writePressureCorrectedVelocityToMomentumSolution(time_arg);
  _velocity_boundary_state_valid = false;
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::computeCellVelocity()
{
  computePressureCorrectedCellVelocity();
}

void
ConservativeDiffuseInterfaceRhieChowMassFlux::addMomentumPredictorExplicitForcing(
    const unsigned int system_i, NumericVector<Number> & rhs) const
{
  mooseAssert(system_i < _momentum_implicit_systems.size(),
              "Momentum component index out of range.");

  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
  {
    if (!hasBlocks(elem_info->subdomain_id()))
      continue;

    const auto & dof_indices = elem_info->dofIndices()[_global_momentum_system_numbers[system_i]];
    if (dof_indices.empty())
      continue;

    const Real force_density =
        reducedPressureMomentumPredictorForceDensity(*elem_info, Moose::currentState())(system_i);
    rhs.add(dof_indices[0], force_density * elem_info->volume() * elem_info->coordFactor());
  }
  rhs.close();
}
