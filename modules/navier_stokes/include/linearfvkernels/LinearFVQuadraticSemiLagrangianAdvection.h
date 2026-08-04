//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "LinearFVElementalKernel.h"

#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

namespace libMesh
{
class PointLocatorBase;
}

/**
 * Advances a cell-centered scalar with a midpoint characteristic trace and a quadratic
 * reconstruction of the old solution at the departure point.
 */
class LinearFVQuadraticSemiLagrangianAdvection : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  LinearFVQuadraticSemiLagrangianAdvection(const InputParameters & params);

  void initialSetup() override;
  void meshChanged() override;

  Real computeMatrixContribution() override;
  Real computeRightHandSideContribution() override;

protected:
  struct Reconstruction
  {
    std::vector<const ElemInfo *> samples;
    std::vector<std::vector<Real>> coefficient_weights;
    Real scale = 0.0;
  };

  void resetGeometry();
  const ElemInfo & containingElement(const Point & point) const;
  const Reconstruction & reconstruction(const ElemInfo & center);
  std::vector<const ElemInfo *> faceNeighbors(const ElemInfo & center) const;
  std::vector<Real> quadraticBasis(const Point & offset, Real scale) const;
  Real reconstructedVariable(const Point & point, const Moose::StateArg & state);
  Real reconstructedFunctor(const Moose::Functor<Real> & functor,
                            const Point & point,
                            const Moose::StateArg & state);
  Point boundedPoint(Point point) const;
  Point departurePoint(const ElemInfo & arrival);

  const unsigned int _dimension;
  const std::array<const Moose::Functor<Real> *, LIBMESH_DIM> _velocity;
  const Real _maximum_courant_number;
  std::array<Real, LIBMESH_DIM> _domain_min;
  std::array<Real, LIBMESH_DIM> _domain_max;
  std::unique_ptr<libMesh::PointLocatorBase> _point_locator;
  std::unordered_map<dof_id_type, Reconstruction> _reconstructions;
};
