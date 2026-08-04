//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "BlockRestrictable.h"
#include "GeneralUserObject.h"
#include "MooseLinearVariableFV.h"
#include "NonADFunctorInterface.h"

#include <array>
#include <unordered_map>
#include <vector>

class LinearSystem;
class ElemInfo;

/**
 * Restores a transported level set to signed distance with an explicit second-order TVD
 * Runge-Kutta method and applies a uniform shift so its regularized-Heaviside material mass equals
 * the conserved reference mass captured before the first correction.
 *
 * A uniform shift leaves grad(phi), normals, and curvature unchanged. Weighting the Heaviside by
 * the material phase density conserves mass while allowing the represented material volume to
 * change when the solid and liquid densities differ.
 */
class MassConservedLevelSetCorrector : public GeneralUserObject,
                                       public NonADFunctorInterface,
                                       public BlockRestrictable
{
public:
  static InputParameters validParams();

  MassConservedLevelSetCorrector(const InputParameters & params);

  void meshChanged() override;
  void initialize() override {}
  void execute() override;
  void finalize() override {}

  Real correct();
  Real correctMassOnly();
  void correctAfterMeshChange();

  const SolverSystemName & systemName() const { return _system_name; }
  const VariableName & variableName() const { return _variable_name; }
  const MooseFunctorName & materialFractionName() const { return _material_fraction_name; }
  const MooseFunctorName & volumetricFaceFluxName() const { return _volumetric_face_flux_name; }
  Real normalizedUpdate() const { return _last_normalized_update; }

protected:
  struct DirectionalStencil
  {
    std::vector<const ElemInfo *> face_neighbors;
    std::vector<Real> face_weights;
    std::vector<Real> reconstruction_weights;
    Real spacing = 0.0;
    bool conforming = false;
  };

  struct CellStencil
  {
    const ElemInfo * elem_info = nullptr;
    std::array<DirectionalStencil, 2 * LIBMESH_DIM> directions;
    std::vector<const ElemInfo *> reconstruction_cells;
    std::array<std::vector<Real>, LIBMESH_DIM> second_derivative_weights;
  };

  using LevelSetValues = std::unordered_map<dof_id_type, Real>;

  Real regularizedHeaviside(Real phi) const;
  Real materialMass(Real shift, const Moose::StateArg & state) const;
  Real massConservingShift(Real target_mass, Real uncorrected_mass) const;
  Real redistance(const Moose::StateArg & state = Moose::currentState());
  void buildRedistanceStencils();
  const CellStencil & stencil(const ElemInfo & elem_info) const;
  const DirectionalStencil &
  direction(const ElemInfo & elem_info, unsigned int axis, bool positive) const;
  Real value(const ElemInfo & elem_info, const LevelSetValues & values) const;
  Real directionalValue(const ElemInfo & elem_info,
                        unsigned int axis,
                        bool positive,
                        const LevelSetValues & values) const;
  Real secondDerivative(const ElemInfo & elem_info,
                        unsigned int axis,
                        const LevelSetValues & values) const;
  Real subcellDistance(const ElemInfo & elem_info,
                       unsigned int axis,
                       bool positive,
                       const LevelSetValues & initial_values) const;
  Real oneSidedDerivative(const ElemInfo & elem_info,
                          unsigned int axis,
                          bool positive,
                          const LevelSetValues & values,
                          const LevelSetValues & initial_values) const;
  Real godunovHamiltonian(const ElemInfo & elem_info,
                          Real initial_value,
                          const LevelSetValues & values,
                          const LevelSetValues & initial_values) const;
  Real localPseudoTimeStep(const ElemInfo & elem_info, const LevelSetValues & initial_values) const;
  void synchronizeLevelSet(LevelSetValues & values, const Moose::StateArg & state);
  void applyShift(Real shift, const Moose::StateArg & state = Moose::currentState());
  void restoreCurrentSolution();

  const SolverSystemName _system_name;
  const VariableName _variable_name;
  const MooseFunctorName _material_fraction_name;
  const MooseFunctorName _volumetric_face_flux_name;
  const Real _interface_width;
  const Real _relative_tolerance;
  const unsigned int _maximum_iterations;
  const unsigned int _redistance_iterations;
  const bool _report;
  Real & _conserved_mass;

  MooseLinearVariableFVReal & _level_set_variable;
  LinearSystem & _system;
  const unsigned int _system_number;
  const unsigned int _variable_number;
  const Moose::Functor<Real> & _material_density;

  std::unordered_map<dof_id_type, CellStencil> _redistance_stencils;
  std::unordered_map<dof_id_type, const ElemInfo *> _redistance_value_cells;
  std::vector<const ElemInfo *> _owned_redistance_cells;
  Real _minimum_grid_spacing = 0.0;
  Real _last_normalized_update = 0.0;
  bool _correcting_after_mesh_change = false;
};
