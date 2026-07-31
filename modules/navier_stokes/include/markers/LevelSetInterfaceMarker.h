//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "Marker.h"
#include "MooseLinearVariableFV.h"

#include <array>
#include <vector>

class LinearSystem;
class MassConservedLevelSetCorrector;

namespace libMesh
{
template <typename>
class NumericVector;
}

/**
 * Refines a cell-width-scaled band around a finite-volume signed-distance interface.
 *
 * A separate, wider coarsening threshold supplies hysteresis. Cells whose centers lie between the
 * two thresholds are left unchanged. A sign change across any face always requests refinement,
 * including at a coarse-fine interface. Before adaptation, limited least-squares gradients are
 * cached in projected system vectors so newly refined children receive conservative second-order
 * cell averages for every available time state.
 */
class LevelSetInterfaceMarker : public Marker
{
public:
  static InputParameters validParams();

  LevelSetInterfaceMarker(const InputParameters & parameters);

  void initialSetup() override;
  void markerSetup() override;
  void meshChanged() override;

protected:
  MarkerValue computeElementMarker() override;

  void initializeProlongationVectors();
  void cacheProlongationGradients();
  void prolongRefinedCells();
  void restrictCoarsenedCells();
  Real regularizedHeaviside(Real phi) const;
  Real conservativeRefinementShift(const std::vector<Real> & child_values,
                                   Real parent_fraction) const;
  Real cellAverageMaterialFraction(Real center_value,
                                   const RealVectorValue & gradient,
                                   const std::array<Real, LIBMESH_DIM> & cell_extent,
                                   unsigned int dimension) const;

  MooseLinearVariableFVReal & _level_set_variable;
  LinearSystem & _level_set_system;
  const unsigned int _system_number;
  const unsigned int _variable_number;
  const Real _refine_band_cells;
  const Real _coarsen_band_cells;
  const bool _conservative_prolongation;
  const Real _conservative_interface_width;
  MassConservedLevelSetCorrector * const _level_set_corrector;
  std::vector<std::vector<libMesh::NumericVector<Number> *>> _prolongation_gradients;
  std::vector<libMesh::NumericVector<Number> *> _projected_material_fractions;
};
