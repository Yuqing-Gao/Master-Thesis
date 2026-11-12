cmake_minimum_required(VERSION 3.20)
project(amg_test LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_executable(amg_test
  src/benchmarks/pattern_experiments.f90
  src/pattern_selection/connection_score.f90
  src/pattern_selection/energy_score.f90
  src/pattern_selection/coupling_strength.f90
  src/pattern_selection/candidate_ranking.f90
  src/matrix_ops.f90
  src/patterns.f90
  src/constraints.f90
  src/minimization.f90
  src/cases/scalar_1d.f90
  src/cases/scalar_2d.f90
  src/main.f90
)

target_compile_options(amg_test PRIVATE
  $<$<CONFIG:Debug>:-g -O0 -Wall -Wextra -fcheck=all -fbacktrace>
  $<$<CONFIG:Release>:-O3>
)