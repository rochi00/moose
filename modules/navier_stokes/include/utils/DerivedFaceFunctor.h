//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FaceCenteredMapFunctor.h"

#include <functional>

/**
 * A face functor whose value follows from other functors rather than from stored data.
 *
 * A quantity that is an exact function of fields the provider already holds does not need a cache
 * of its own: a cache of it can only ever agree or disagree with its definition, and keeping the
 * two in step means remembering to refresh it whenever an input moves. Expressing it as a function
 * instead removes that obligation, and with it the possibility of reading a value that was correct
 * for an earlier state.
 *
 * Face evaluation is the only thing this redefines. Cell reconstruction, the face argument
 * overloads and the block restriction are all inherited, so a derived quantity behaves exactly
 * like a stored one everywhere it is consumed.
 */
template <typename T>
class DerivedFaceFunctor
  : public FaceCenteredMapFunctor<T, std::unordered_map<dof_id_type, T>>
{
public:
  using Base = FaceCenteredMapFunctor<T, std::unordered_map<dof_id_type, T>>;
  using typename Base::ValueType;

  /// @param face_value the definition of this quantity on a face
  DerivedFaceFunctor(const MooseMesh & mesh,
                     const std::string & name,
                     std::function<T(const FaceInfo &)> face_value)
    : Base(mesh, name), _face_value(std::move(face_value))
  {
  }

  DerivedFaceFunctor(const MooseMesh & mesh,
                     const std::set<SubdomainID> & sub_ids,
                     const std::string & name,
                     std::function<T(const FaceInfo &)> face_value)
    : Base(mesh, sub_ids, name), _face_value(std::move(face_value))
  {
  }

  ValueType evaluate(const FaceInfo * const fi) const override
  {
    mooseAssert(fi, "A derived face functor needs a face.");
    mooseAssert(this->empty(),
                "A derived face functor is defined by its function, so nothing should have been "
                "stored into it.");
    return _face_value(*fi);
  }

private:
  const std::function<T(const FaceInfo &)> _face_value;
};
