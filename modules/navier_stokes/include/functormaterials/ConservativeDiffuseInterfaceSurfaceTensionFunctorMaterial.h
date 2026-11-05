//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"

#include <set>

class MooseMesh;
class ElemInfo;

/**
 * Publishes the continuum surface force of a diffuse interface: the interface curvature, and the
 * surface tension force density sigma * kappa * grad(alpha) that a momentum equation consumes.
 *
 * The curvature is the divergence of the interface unit normal, taken over the faces of a cell,
 * which is why this is a least-squares reconstruction rather than a functor gradient: the normal
 * has to be assembled per face before being differenced. On the faces of the contact-angle
 * boundaries the normal is the one a static contact angle prescribes, which is how the contact
 * angle enters the curvature of the wall cells.
 */
class ConservativeDiffuseInterfaceSurfaceTensionFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  ConservativeDiffuseInterfaceSurfaceTensionFunctorMaterial(const InputParameters & parameters);

protected:
  /// The interface curvature of one cell, as the divergence of the unit normal over its faces.
  Real cellCurvature(const ElemInfo & elem_info, const Moose::StateArg & time_arg) const;

  /// The surface tension force density of one cell.
  RealVectorValue cellForceDensity(const ElemInfo & elem_info,
                                   const Moose::StateArg & time_arg) const;

  /// Curvature at an argument that names a single cell.
  template <typename SpaceArg>
  Real curvature(const SpaceArg & r, const Moose::StateArg & t) const
  {
    return cellCurvature(_moose_mesh.elemInfo(r.elem->id()), t);
  }

  /// Curvature on a face, interpolated from the cells it couples.
  Real curvature(const Moose::FaceArg & r, const Moose::StateArg & t) const;

  /// Both quantities are cell based, so there is no nodal evaluation of them.
  Real curvature(const Moose::NodeArg &, const Moose::StateArg &) const;

  /// Force density at an argument that names a single cell.
  template <typename SpaceArg>
  RealVectorValue forceDensity(const SpaceArg & r, const Moose::StateArg & t) const
  {
    return cellForceDensity(_moose_mesh.elemInfo(r.elem->id()), t);
  }

  /// Force density on a face, interpolated from the cells it couples.
  RealVectorValue forceDensity(const Moose::FaceArg & r, const Moose::StateArg & t) const;

  RealVectorValue forceDensity(const Moose::NodeArg &, const Moose::StateArg &) const;

  /// Whether a face couples two cells this material is defined on.
  bool isInternalFace(const FaceInfo & fi) const;
  const ElemInfo & boundaryElemInfo(const FaceInfo & fi) const;

  /// Whether a boundary face lies on a contact-angle boundary
  bool onContactAngleBoundary(const FaceInfo & fi) const;
  /// The interface normal on a wall face with the static contact angle, from the cell's normal
  /// (pointing into the liquid, of the magnitude a vanishing gradient regularizes to) and the
  /// unit wall normal pointing out of the domain
  RealVectorValue contactAngleNormal(const RealVectorValue & cell_normal,
                                     const RealVectorValue & wall_normal) const;

  /// The volume fraction whose interface this describes
  const Moose::Functor<Real> & _volume_fraction;

  /// The surface tension coefficient
  const Moose::Functor<Real> & _surface_tension_coefficient;

  /// Added to |grad(alpha)| when forming the unit normal, so a vanishing gradient is finite
  const Real _normal_regularization;

  /// Cosine and sine of the static contact angle on the contact-angle boundaries
  const Real _cos_contact_angle;
  const Real _sin_contact_angle;
  std::set<BoundaryID> _contact_angle_boundaries;

  const MooseMesh & _moose_mesh;
  const unsigned int _dim;
};
