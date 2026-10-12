cmake_minimum_required(VERSION 3.20)
project(fv_layout_bench CXX)

# Uses the Kokkos installed by PETSc in the MOOSE conda environment (or any Kokkos install):
#   cmake -S . -B build -DKokkos_ROOT=$CONDA_PREFIX -DCMAKE_BUILD_TYPE=Release
find_package(Kokkos REQUIRED)

add_executable(fv_layout_bench fv_layout_bench.cpp)
target_link_libraries(fv_layout_bench Kokkos::kokkos)
target_compile_features(fv_layout_bench PRIVATE cxx_std_17)
