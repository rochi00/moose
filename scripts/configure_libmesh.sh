#!/usr/bin/env bash
#* This file is part of the MOOSE framework
#* https://mooseframework.inl.gov
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

# Configure libMesh with the default MOOSE configuration options
#
# Separated so that it can be used across all libMesh build scripts:
# - scripts/update_and_rebuild_libmesh.sh
# - conda/libmesh/build.sh
function configure_libmesh()
{
  if [ -z "$SRC_DIR" ]; then
    echo "SRC_DIR is not set for configure_libmesh"
    exit 1
  fi

  if [ ! -d "$SRC_DIR" ]; then
    echo "$SRC_DIR=SRC_DIR does not exist"
    exit 1
  fi

  if [ -z "$LIBMESH_DIR" ]; then
    echo "$LIBMESH_DIR is not set for configure_libmesh"
    exit 1
  fi

  # Preserves capability in update_and_rebuild_libmesh.sh, but this is set in
  # conda/libmesh/build.sh. If not, conda considers it an "unbound variable"
  if [[ -n "$INSTALL_BINARY" ]]; then
    echo "INFO: INSTALL_BINARY set"
  else
    export INSTALL_BINARY="${SRC_DIR}/build-aux/install-sh -C"
  fi

  # If METHODS is not set in update_and_rebuild_libmesh.sh, set a default value.
  export METHODS=${METHODS:="opt oprof devel dbg"}

  EXTRA_ARGS=()
  # libtirpc has changed paths from a previous default searched location
  if [[ $(uname) == Linux ]] && [[ -d ${CONDA_PREFIX}/include/tirpc ]]; then
    EXTRA_ARGS+=("--with-xdr-include=${CONDA_PREFIX}/include/tirpc")
  fi
  # Enable native Kokkos FE math headers if Kokkos is available via PETSc.
  # Pass the same compiler/flags that MOOSE's kokkos.mk uses so that
  # libMesh's .K test files compile identically.
  if [[ -n "$PETSC_DIR" ]] && [[ -f "${PETSC_DIR}/include/Kokkos_Core.hpp" ]]; then
    EXTRA_ARGS+=("--with-kokkos=${PETSC_DIR}")

    KOKKOS_CFG="${PETSC_DIR}/include/KokkosCore_config.h"
    PETSC_VARS="${PETSC_DIR}/lib/petsc/conf/petscvariables"
    _petsc_makevar()
    {
      local var_name="$1"
      [[ -r "$PETSC_VARS" ]] || return 0
      sed -n "s/^${var_name} = //p" "$PETSC_VARS" | head -n 1
    }
    _kokkos_openmp=""
    _kokkos_has_cuda=0
    _kokkos_has_hip=0
    _kokkos_has_sycl=0
    if [[ -r "$KOKKOS_CFG" ]]; then
      grep -q '^#define KOKKOS_ENABLE_OPENMP' "$KOKKOS_CFG" && _kokkos_openmp="-fopenmp"
      grep -q '^#define KOKKOS_ENABLE_CUDA' "$KOKKOS_CFG" && _kokkos_has_cuda=1
      grep -q '^#define KOKKOS_ENABLE_HIP' "$KOKKOS_CFG" && _kokkos_has_hip=1
      grep -q '^#define KOKKOS_ENABLE_SYCL' "$KOKKOS_CFG" && _kokkos_has_sycl=1
    fi

    # Detect backend and arch from PETSc's petscconf.h (mirrors kokkos.mk logic)
    _petsc_conf="${PETSC_DIR}/include/petscconf.h"
    if [[ -n "${PETSC_ARCH:-}" ]] && [[ -r "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h" ]]; then
      _petsc_conf="${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h"
    fi
    _petsc_have_cuda=$(sed -n 's/#define PETSC_HAVE_CUDA //p' "$_petsc_conf" 2>/dev/null)
    _petsc_have_hip=$(sed -n 's/#define PETSC_HAVE_HIP //p' "$_petsc_conf" 2>/dev/null)
    _petsc_have_sycl=$(sed -n 's/#define PETSC_HAVE_SYCL //p' "$_petsc_conf" 2>/dev/null)
    _use_cuda=0
    _use_hip=0
    _use_sycl=0

    [[ "$_petsc_have_cuda" == "1" || "$_kokkos_has_cuda" == "1" ]] && _use_cuda=1
    [[ "$_petsc_have_hip" == "1" || "$_kokkos_has_hip" == "1" ]] && _use_hip=1
    [[ "$_petsc_have_sycl" == "1" || "$_kokkos_has_sycl" == "1" ]] && _use_sycl=1

    if [[ "$_use_cuda" == "1" ]] && command -v nvcc &>/dev/null; then
      # For CUDA, let libMesh choose nvcc_wrapper itself. Exporting raw
      # KOKKOS_CXX=nvcc here overrides libMesh's safer CUDA toolchain logic
      # and causes ordinary host flags to hit nvcc project-wide.
      unset KOKKOS_CXX KOKKOS_CXXFLAGS KOKKOS_LDFLAGS

    elif [[ "$_use_hip" == "1" ]] && command -v hipcc &>/dev/null; then
      export KOKKOS_CXX="$(command -v hipcc)"
      export KOKKOS_CXXFLAGS="${_kokkos_openmp}"
      export KOKKOS_LDFLAGS="-L${PETSC_DIR}/lib"

    elif [[ "$_use_sycl" == "1" ]] && command -v icpx &>/dev/null; then
      export KOKKOS_CXX="$(command -v icpx)"
      export KOKKOS_CXXFLAGS="-fsycl ${_kokkos_openmp}"
      export KOKKOS_LDFLAGS="-L${PETSC_DIR}/lib"

    else
      # CPU-only (OpenMP or serial)
      export KOKKOS_CXX="${CXX}"
      export KOKKOS_CXXFLAGS="${_kokkos_openmp} -x c++"
      export KOKKOS_LDFLAGS="-L${PETSC_DIR}/lib"
      EXTRA_ARGS+=("--with-kokkos-backend=openmp")
    fi

    if [[ -z "${KOKKOS_CPPFLAGS:-}" ]]; then
      export KOKKOS_CPPFLAGS="-DLIBMESH_KOKKOS_COMPILATION -I${PETSC_DIR}/include"
    fi

    if [[ -z "${KOKKOS_LIBS:-}" ]]; then
      _petsc_kokkos_libs="$(_petsc_makevar KOKKOS_LIB)"
      _petsc_backend_libs=""
      _petsc_system_libs=""

      if [[ "$_use_cuda" == "1" ]]; then
        _petsc_backend_libs="$(_petsc_makevar CUDA_LIB)"
      elif [[ "$_use_hip" == "1" ]]; then
        _petsc_backend_libs="$(_petsc_makevar HIP_LIB)"
      elif [[ "$_use_sycl" == "1" ]]; then
        _petsc_backend_libs="$(_petsc_makevar SYCL_LIB)"
      fi

      _petsc_external_libs="$(_petsc_makevar PETSC_EXTERNAL_LIB_BASIC)"
      case " $_petsc_external_libs " in
        *" -ldl "*) _petsc_system_libs="$_petsc_system_libs -ldl" ;;
      esac
      case " $_petsc_external_libs " in
        *" -lpthread "*) _petsc_system_libs="$_petsc_system_libs -lpthread" ;;
      esac

      if [[ -n "$_petsc_kokkos_libs" ]]; then
        export KOKKOS_LIBS="$_petsc_kokkos_libs $_petsc_backend_libs $_petsc_system_libs"
      else
        export KOKKOS_LIBS="-lkokkoscore $_petsc_system_libs"
      fi
    fi
  fi

  # Allow unbound variable for when EXTRA_ARGS is empty
  set +u

  cd "${SRC_DIR}/build" || exit 1
  # shellcheck disable=SC2086  # we want wordsplitting
  # shellcheck disable=SC2048  # we want to not handle whitspaces
  ../configure --enable-silent-rules \
               --enable-unique-id \
               --disable-warnings \
               --with-thread-model=openmp \
               --disable-maintainer-mode \
               --enable-hdf5 \
               --enable-petsc-hypre-required \
               --enable-metaphysicl-required \
               --enable-xdr-required \
               --with-cxx-std-min=2014 \
               --without-gdb-command \
               --with-methods="${METHODS}" \
               --prefix="${LIBMESH_DIR}" \
               --with-future-timpi-dir="${LIBMESH_DIR}" \
               INSTALL="${INSTALL_BINARY}" \
               "${EXTRA_ARGS[@]}" \
               $*
  local RETURN_CODE=$?

  # Restore unbound variable checks
  set -u

  return $RETURN_CODE
}
