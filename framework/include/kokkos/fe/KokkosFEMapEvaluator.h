//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosArray.h"
#include "libmesh/gpu/kokkos_fe_types.h"
#include "libmesh/gpu/kokkos_fe_lagrange_1d.h"
#include "libmesh/gpu/kokkos_fe_lagrange_2d.h"
#include "libmesh/gpu/kokkos_fe_lagrange_3d.h"

namespace Moose::Kokkos
{
using namespace libMesh::Kokkos;

/**
 * CPU-side evaluator for geometric map arrays.
 *
 * Provides static helpers that fill the four geometric-map Array2D objects in
 * KokkosAssembly::initShape() without calling libMesh FEBase::reinit():
 *
 *   _map_phi            parent LAGRANGE shape functions at volume quad points
 *   _map_grad_phi       their reference-domain gradients
 *   _map_psi_face       side-element LAGRANGE shape functions at face quad points
 *   _map_grad_psi_face  their reference-domain gradients (dpsidxi, dpsideta)
 *
 * Each function both creates and fills its destination array.  All four must be
 * called from CPU code only (they are compiled only when MOOSE_KOKKOS_SCOPE is
 * defined, i.e. from .K translation units).
 *
 * Coordinate conventions match libMesh:
 *   - qpts  are in the parent element's reference coordinate system
 *   - face_qpts  are in the side element's reference coordinate system
 *     (as stored in _q_points_face, which is initialised by qrule_face->init(*side_elem))
 */
class FEMapEvaluator
{
public:
#ifdef MOOSE_KOKKOS_SCOPE

  /**
   * Fill dest = _map_phi(sid, etid) for parent topology @p topo.
   *
   * dest(node, qp) = phi_node(qpts[qp])
   *
   * @param dest   Array2D<Real> to create and fill (node × qp)
   * @param topo   Parent element topology
   * @param qpts   Volume quad points in the parent reference coordinate system
   */
  static inline void
  fillMapPhi(Array2D<Real> & dest, libMesh::ElemType topo, const Array<Real3> & qpts)
  {
    const unsigned int n_qp = qpts.size();
    switch (topo)
    {
      case libMesh::EDGE2:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE2>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE2>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::EDGE3:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE3>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE3>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::TRI3:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::TRI6:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::QUAD4:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::QUAD8:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::QUAD9:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::TET4:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TET4>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::TET4>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::TET10:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TET10>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::TET10>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::HEX8:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::HEX8>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::HEX8>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::HEX20:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::HEX20>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::HEX20>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::HEX27:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::HEX27>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) =
                FEEvaluator<libMesh::LAGRANGE, libMesh::HEX27>::shape(i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      default:
        mooseError("FEMapEvaluator::fillMapPhi: unsupported element topology");
    }
  }

  /**
   * Fill dest = _map_grad_phi(sid, etid) for parent topology @p topo.
   *
   * dest(node, qp) = grad_phi_node(qpts[qp])  (reference-domain gradient)
   *
   * @param dest   Array2D<Real3> to create and fill (node × qp)
   * @param topo   Parent element topology
   * @param qpts   Volume quad points in the parent reference coordinate system
   */
  static inline void
  fillMapGradPhi(Array2D<Real3> & dest, libMesh::ElemType topo, const Array<Real3> & qpts)
  {
    const unsigned int n_qp = qpts.size();
    switch (topo)
    {
      case libMesh::EDGE2:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE2>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE2>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::EDGE3:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE3>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE3>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::TRI3:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::TRI6:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::QUAD4:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::QUAD8:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::QUAD9:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::TET4:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TET4>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::TET4>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::TET10:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TET10>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::TET10>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::HEX8:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::HEX8>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::HEX8>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::HEX20:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::HEX20>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::HEX20>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      case libMesh::HEX27:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::HEX27>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::HEX27>::grad_shape(
                i, qpts[qp](0), qpts[qp](1), qpts[qp](2));
        break;
      }
      default:
        mooseError("FEMapEvaluator::fillMapGradPhi: unsupported element topology");
    }
  }

  /**
   * Fill dest = _map_psi_face(sid, etid)(side) for side topology @p side_topo.
   *
   * dest(node, qp) = psi_node(face_qpts[qp])
   *
   * face_qpts are the face quadrature points expressed in the side element's own
   * reference coordinate system (as stored in _q_points_face, which is populated
   * by qrule_face->init(*elem->side_ptr(side))).
   *
   * @param dest       Array2D<Real> to create and fill (side_node × face_qp)
   * @param side_topo  Side element topology (from getSideTopology(parent_topo))
   * @param face_qpts  Face quad points in the side element's reference coordinate system
   */
  static inline void
  fillMapPsiFace(Array2D<Real> & dest,
                 libMesh::ElemType side_topo,
                 const Array<Real3> & face_qpts)
  {
    const unsigned int n_qp = face_qpts.size();
    switch (side_topo)
    {
      case libMesh::EDGE2:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE2>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE2>::shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
        break;
      }
      case libMesh::EDGE3:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE3>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE3>::shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
        break;
      }
      case libMesh::TRI3:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>::shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
        break;
      }
      case libMesh::TRI6:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
        break;
      }
      case libMesh::QUAD4:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
        break;
      }
      case libMesh::QUAD8:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>::shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
        break;
      }
      case libMesh::QUAD9:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
            dest(i, qp) = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>::shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
        break;
      }
      default:
        mooseError("FEMapEvaluator::fillMapPsiFace: unsupported side topology");
    }
  }

  /**
   * Fill dest = _map_grad_psi_face(sid, etid)(side) for side topology @p side_topo.
   *
   * Matches what libMesh stores in get_dpsidxi / get_dpsideta:
   *   dest(node, qp)(0) = dpsidxi   (set when parent_dim > 1)
   *   dest(node, qp)(1) = dpsideta  (set when parent_dim > 2)
   *   dest(node, qp)(2) = 0         (never used)
   *
   * face_qpts are in the side element's reference coordinate system.
   *
   * @param dest        Array2D<Real3> to create and fill (side_node × face_qp)
   * @param side_topo   Side element topology (from getSideTopology(parent_topo))
   * @param face_qpts   Face quad points in the side element's reference coordinate system
   * @param parent_dim  Spatial dimension of the parent element
   */
  static inline void
  fillMapGradPsiFace(Array2D<Real3> & dest,
                     libMesh::ElemType side_topo,
                     const Array<Real3> & face_qpts,
                     unsigned int parent_dim)
  {
    const unsigned int n_qp = face_qpts.size();
    switch (side_topo)
    {
      case libMesh::EDGE2:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE2>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
          {
            auto g = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE2>::grad_shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
            if (parent_dim > 1)
              dest(i, qp)(0) = g(0);
            if (parent_dim > 2)
              dest(i, qp)(1) = g(1);
          }
        break;
      }
      case libMesh::EDGE3:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE3>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
          {
            auto g = FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE3>::grad_shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
            if (parent_dim > 1)
              dest(i, qp)(0) = g(0);
            if (parent_dim > 2)
              dest(i, qp)(1) = g(1);
          }
        break;
      }
      case libMesh::TRI3:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
          {
            auto g = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>::grad_shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
            if (parent_dim > 1)
              dest(i, qp)(0) = g(0);
            if (parent_dim > 2)
              dest(i, qp)(1) = g(1);
          }
        break;
      }
      case libMesh::TRI6:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
          {
            auto g = FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::grad_shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
            if (parent_dim > 1)
              dest(i, qp)(0) = g(0);
            if (parent_dim > 2)
              dest(i, qp)(1) = g(1);
          }
        break;
      }
      case libMesh::QUAD4:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
          {
            auto g = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::grad_shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
            if (parent_dim > 1)
              dest(i, qp)(0) = g(0);
            if (parent_dim > 2)
              dest(i, qp)(1) = g(1);
          }
        break;
      }
      case libMesh::QUAD8:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
          {
            auto g = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>::grad_shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
            if (parent_dim > 1)
              dest(i, qp)(0) = g(0);
            if (parent_dim > 2)
              dest(i, qp)(1) = g(1);
          }
        break;
      }
      case libMesh::QUAD9:
      {
        constexpr unsigned int n = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>::n_dofs();
        dest.create(n, n_qp);
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int qp = 0; qp < n_qp; ++qp)
          {
            auto g = FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>::grad_shape(
                i, face_qpts[qp](0), face_qpts[qp](1), face_qpts[qp](2));
            if (parent_dim > 1)
              dest(i, qp)(0) = g(0);
            if (parent_dim > 2)
              dest(i, qp)(1) = g(1);
          }
        break;
      }
      default:
        mooseError("FEMapEvaluator::fillMapGradPsiFace: unsupported side topology");
    }
  }

#endif // MOOSE_KOKKOS_SCOPE
};

} // namespace Moose::Kokkos
