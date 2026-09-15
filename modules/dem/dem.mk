# ArborX (header-only) for the bounding-volume-hierarchy broad phase; optional, present when the
# submodule is checked out. ArborX 2.x needs C++20, so the module's Kokkos translation units,
# the only ones that include it, are compiled as C++20 (the flag follows libMesh's -std=c++17)
ARBORX_DIR ?= $(MOOSE_DIR)/modules/dem/contrib/arborx
ifneq ($(wildcard $(ARBORX_DIR)/src/spatial/ArborX_LinearBVH.hpp),)
  ADDITIONAL_CPPFLAGS        += -DDEM_HAVE_ARBORX
  ADDITIONAL_KOKKOS_CPPFLAGS += -std=c++20
  ADDITIONAL_INCLUDES        += -isystem $(ARBORX_DIR)/src -isystem $(ARBORX_DIR)/src/geometry \
                                -isystem $(ARBORX_DIR)/src/spatial \
                                -isystem $(MOOSE_DIR)/modules/dem/contrib/arborx_config
endif
