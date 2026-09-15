// Stand-in for the ArborX_Config.hpp that ArborX's CMake build would generate; the DEM module
// uses the header-only library without CMake, and the DEM broad phase needs no optional feature
// (MPI, rocThrust, oneDPL)
#ifndef ARBORX_CONFIG_HPP
#define ARBORX_CONFIG_HPP

#endif
