//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeLinearFVLimitedGradientThread.h"

#include "GradientLimiterType.h"
#include "SystemBase.h"
#include "PetscVectorReader.h"
#include "FEProblemBase.h"
#include "FVUtils.h"

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_object.h"

#include <algorithm>
#include <cmath>
#include <limits>

ComputeLinearFVLimitedGradientThread::ComputeLinearFVLimitedGradientThread(
    FEProblemBase & fe_problem,
    SystemBase & system,
    const std::vector<std::unique_ptr<NumericVector<Number>>> & raw_gradient,
    std::vector<std::unique_ptr<NumericVector<Number>>> & temporary_limited_gradient,
    const Moose::FV::GradientLimiterType limiter_type,
    const std::unordered_set<unsigned int> & requested_variables)
  : _fe_problem(fe_problem),
    _dim(_fe_problem.mesh().dimension()),
    _system(system),
    _libmesh_system(system.system()),
    _system_number(_libmesh_system.number()),
    _raw_gradient(raw_gradient),
    _limiter_type(limiter_type),
    _requested_variables(requested_variables),
    _temporary_limited_gradient(temporary_limited_gradient)
{
}

ComputeLinearFVLimitedGradientThread::ComputeLinearFVLimitedGradientThread(
    ComputeLinearFVLimitedGradientThread & x, Threads::split /*split*/)
  : _fe_problem(x._fe_problem),
    _dim(x._dim),
    _system(x._system),
    _libmesh_system(x._libmesh_system),
    _system_number(x._system_number),
    _raw_gradient(x._raw_gradient),
    _limiter_type(x._limiter_type),
    _requested_variables(x._requested_variables),
    _temporary_limited_gradient(x._temporary_limited_gradient)
{
}

void
ComputeLinearFVLimitedGradientThread::operator()(const ElemInfoRange & range)
{
  ParallelUniqueId puid;
  _tid = puid.id;

  const auto & raw_grad_container = _raw_gradient;

  if (_limiter_type != Moose::FV::GradientLimiterType::Venkatakrishnan)
    mooseError("ComputeLinearFVLimitedGradientThread currently supports only the Venkatakrishnan "
               "limiter.");

  unsigned int size = 0;

  for (const auto & variable : _system.getVariables(_tid))
  {
    _current_var = dynamic_cast<MooseLinearVariableFV<Real> *>(variable);
    if (!_current_var)
      continue;

    if (!_current_var->needsGradientVectorStorage())
      continue;

    if (!_requested_variables.count(_current_var->number()))
      continue;

    if (!size)
      size = range.size();

    std::vector<std::vector<Real>> temporary_values(_temporary_limited_gradient.size(),
                                                    std::vector<Real>(size, 0.0));
    std::vector<dof_id_type> dof_indices(size, 0);

    PetscVectorReader solution_reader(*_libmesh_system.current_local_solution);
    std::vector<PetscVectorReader> grad_reader;
    grad_reader.reserve(raw_grad_container.size());
    for (const auto dim_index : index_range(raw_grad_container))
      grad_reader.emplace_back(*raw_grad_container[dim_index]);

    mooseAssert(raw_grad_container.size() >= _dim,
                "Raw gradient container has fewer components than mesh dimension.");
    mooseAssert(_temporary_limited_gradient.size() >= _dim,
                "Limited gradient container has fewer components than mesh dimension.");

    auto elem_iterator = range.begin();
    for (const auto elem_i : make_range(size))
    {
      const auto & elem_info = *elem_iterator;
      elem_iterator++;

      if (!_current_var->hasBlocks(elem_info->subdomain_id()))
        continue;

      const dof_id_type dof = elem_info->dofIndices()[_system_number][_current_var->number()];
      if (dof == libMesh::DofObject::invalid_id)
        continue;

      dof_indices[elem_i] = dof;

      const Elem * const elem = elem_info->elem();
      const Real phi_elem = solution_reader(dof);
      const Real scale = elem->hmin();
      Real max_value = phi_elem;
      Real min_value = phi_elem;
      DenseMatrix<Real> normal_matrix(_dim, _dim);
      DenseVector<Real> right_hand_side(_dim);
      std::vector<std::pair<Point, Real>> neighbor_samples;
      std::vector<Point> internal_face_points;

      auto gather_neighbor = [this,
                              &solution_reader,
                              &elem_info,
                              phi_elem,
                              scale,
                              &max_value,
                              &min_value,
                              &normal_matrix,
                              &right_hand_side,
                              &neighbor_samples,
                              &internal_face_points](const Elem &,
                                                     const Elem * const neighbor,
                                                     const FaceInfo * const face_info,
                                                     const Point & surface_vector,
                                                     const Real,
                                                     const bool)
      {
        if (!neighbor || !_current_var->hasBlocks(neighbor->subdomain_id()))
          return;

        const auto & neighbor_info = _fe_problem.mesh().elemInfo(neighbor->id());
        const dof_id_type neighbor_dof =
            neighbor_info.dofIndices()[_system_number][_current_var->number()];
        if (neighbor_dof == libMesh::DofObject::invalid_id)
          return;

        const Real phi_neighbor = solution_reader(neighbor_dof);
        max_value = std::max(max_value, phi_neighbor);
        min_value = std::min(min_value, phi_neighbor);

        const Point normalized_offset = (neighbor_info.centroid() - elem_info->centroid()) / scale;
        const Real offset_norm_sq = normalized_offset.norm_sq();
        if (offset_norm_sq <= std::numeric_limits<Real>::epsilon())
          return;

        // Face measure prevents the multiple fine children on one split coarse face from
        // receiving more aggregate influence merely because that face has been subdivided.
        const Real weight = surface_vector.norm() / offset_norm_sq;
        const Real value_delta = phi_neighbor - phi_elem;
        for (const auto row : make_range(_dim))
        {
          right_hand_side(row) += weight * normalized_offset(row) * value_delta;
          for (const auto column : make_range(_dim))
            normal_matrix(row, column) +=
                weight * normalized_offset(row) * normalized_offset(column);
        }

        neighbor_samples.emplace_back(normalized_offset, value_delta);
        internal_face_points.push_back(face_info->faceCentroid());
      };

      const auto coordinate_system = _fe_problem.mesh().getCoordSystem(elem_info->subdomain_id());
      Moose::FV::loopOverElemFaceInfo(*elem,
                                      _fe_problem.mesh(),
                                      gather_neighbor,
                                      coordinate_system,
                                      _fe_problem.mesh().getAxisymmetricRadialCoord());

      // Keep the preexisting Green-Gauss result only as a fallback for a rank-deficient stencil.
      VectorValue<Real> raw_grad;
      raw_grad.zero();
      for (const auto dim_index : make_range(_dim))
        raw_grad(dim_index) = grad_reader[dim_index](dof);

      if (neighbor_samples.size() >= _dim)
      {
        DenseVector<Real> scaled_gradient(_dim);
        normal_matrix.svd_solve(right_hand_side, scaled_gradient);
        for (const auto dim_index : make_range(_dim))
          raw_grad(dim_index) = scaled_gradient(dim_index) / scale;
      }

      // If the stencil is constant (or nearly constant), don't attempt to limit.
      if (std::abs(max_value - min_value) < 1e-14)
      {
        for (const auto dim_index : make_range(_dim))
          temporary_values[dim_index][elem_i] = raw_grad(dim_index);
        continue;
      }

      Real alpha = 1.0;
      const Point & elem_centroid = elem_info->centroid();

      // An affine field is represented exactly by the least-squares polynomial. Do not let the
      // nonlinear limiter destroy that consistency on unequal coarse-fine center spacings.
      Real maximum_fit_residual = 0.0;
      Real value_scale = std::max(1.0, std::abs(phi_elem));
      for (const auto & [normalized_offset, value_delta] : neighbor_samples)
      {
        maximum_fit_residual = std::max(
            maximum_fit_residual, std::abs(value_delta - raw_grad * (scale * normalized_offset)));
        value_scale = std::max(value_scale, std::abs(phi_elem + value_delta));
      }

      const bool affine_stencil = maximum_fit_residual <= 1e-12 * value_scale;

      if (!affine_stencil)
        for (const auto & face_point : internal_face_points)
        {
          const Real delta_face = raw_grad * (face_point - elem_centroid);
          const Real admissible_delta =
              delta_face >= 0.0 ? max_value - phi_elem : phi_elem - min_value;
          if (std::abs(delta_face) <= admissible_delta * (1.0 + 1e-12))
            continue;

          if (admissible_delta <= 0.0)
          {
            alpha = 0.0;
            continue;
          }

          const Real ratio = std::abs(delta_face) / admissible_delta;
          const Real beta = (2.0 * ratio + 1.0) / (ratio * (2.0 * ratio + 1.0) + 1.0);
          alpha = std::min(alpha, beta);
        }

      const VectorValue<Real> limited_grad = alpha * raw_grad;
      for (const auto dim_index : make_range(_dim))
        temporary_values[dim_index][elem_i] = limited_grad(dim_index);
    }

    for (const auto dim_index : make_range(_dim))
    {
      _temporary_limited_gradient[dim_index]->add_vector(temporary_values[dim_index].data(),
                                                         dof_indices);
    }
  }
}

void
ComputeLinearFVLimitedGradientThread::join(const ComputeLinearFVLimitedGradientThread & /*y*/)
{
}
