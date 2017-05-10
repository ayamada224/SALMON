# recommended version : 3.0.2
cmake_minimum_required(VERSION 2.8)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/misc.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/setup.cmake)


### Project settings
project(SALMON Fortran C)

if (${CMAKE_CROSSCOMPILING})
  #
  # (For supercomputer systems) cross compile mode, use toolchain file.
  #
  message(STATUS "Target platform ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
else ()
  # (For PC clusters and PC) compiler environment detects automatically.
  message(STATUS "Target platform is Native (${CMAKE_HOST_SYSTEM})")

  set(TARGET_SUFFIX ".cpu")

  check_mpi_compiler(${CMAKE_Fortran_COMPILER} IS_MPI_COMPILER)
  if (${IS_MPI_COMPILER})
    set(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})
  endif ()

  check_mpi_compiler(${CMAKE_C_COMPILER} IS_MPI_COMPILER)
  if (${IS_MPI_COMPILER})
    set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
  endif ()

  find_package(MPI REQUIRED)

  if (NOT DEFINED MPI_Fortran_FOUND)
    message(FATAL_ERROR "MPI Fortran compilers not found.")
  endif()

  if (NOT DEFINED MPI_C_FOUND)
    message(FATAL_ERROR "MPI C compilers not found.")
  endif()

  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

  set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(Fortran_FLAGS_General       "-cpp ${MPI_Fortran_COMPILE_FLAGS}")

  set(CMAKE_C_FLAGS_DEBUG         "-O2 -g")
  set(CMAKE_C_FLAGS_RELEASE       "-O3")
  set(C_FLAGS_General             "${MPI_C_COMPILE_FLAGS}")

  find_package(OpenMP    REQUIRED)
  find_package(ScaLAPACK REQUIRED)

  set(OPENMP_FLAGS    ${OpenMP_C_FLAGS})
  set(ScaLAPACK_FLAGS ${ScaLAPACK_LINKER_FLAGS} ${ScaLAPACK_LIBRARIES})
endif ()

message(STATUS "user-defined FFLAGS = ${CMAKE_Fortran_FLAGS}")
message(STATUS "user-defined CFLAGS = ${CMAKE_C_FLAGS}")


### Set optimize options
option_set(USE_TLOG                       "Use TLOG library"                                         OFF)
option_set(USE_PAPI                       "Use Performance API"                                      OFF)
option_set(USE_NVTX                       "Use NVIDIA Tools Extension Library"                       OFF)
option_set(OPT_STENCIL                    "Optimize stencil code"                                    ON)
option_set(STENCIL_WITH_C                 "Use C code at stencil"                                    OFF)
option_set(DOMAIN_IS_POW2                 "Stencil domain is power of two"                           OFF)
option_set(ENABLE_CURRENT_PREPROCESSING   "Enable preprocessing of the current computation in RT"    OFF)
option_set(ENABLE_ARRAY_PADDING           "Enable array padding for the stencil"                     ON)
option_set(ENABLE_EXPLICIT_VEC            "Enable explicit vectorization for the stencil"            OFF)
option_set(ENABLE_LOOP_BLOCKING           "Enable loop blocking at the stencil"                      OFF)
option_set(ENABLE_SWPREFETCH              "Enable software prefetch in the explicit vec of stencil"  OFF)
option_set(ENABLE_REDUCE_FOR_MANYCORE     "Enable optimized reduction code for many-core processor"  OFF)
option_set(ENABLE_LARGE_BLOCKING          "Enable large blocking for accelerator"                    OFF)
option_set(USE_OLD_PROPAGATOR             "Use old propagator for comparing to past results"         OFF)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}")
set(TARGET_NAME            "salmon${TARGET_SUFFIX}")

set(CMAKE_Fortran_FLAGS    "${ARCH} ${OPENMP_FLAGS} ${Fortran_FLAGS_General} ${ADDITIONAL_OPTIMIZE_FLAGS} ${CMAKE_Fortran_FLAGS}")
set(CMAKE_C_FLAGS          "${ARCH} ${OPENMP_FLAGS} ${C_FLAGS_General} ${ADDITIONAL_OPTIMIZE_FLAGS} ${CMAKE_C_FLAGS}")
set(EXTERNAL_LIBS          "${ScaLAPACK_FLAGS}")

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

### Performance API from ICL, Univ. TENNESSEE
if (USE_PAPI)
  add_definitions(-DARTED_USE_PAPI)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} papi)
endif ()

### NVIDIA Tools Extension Library
if (USE_NVTX)
  find_package(CUDA REQUIRED)
  add_definitions(-DARTED_USE_NVTX)
  link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} nvToolsExt)
endif ()

add_definitions_if(ENABLE_CURRENT_PREPROCESSING -DARTED_CURRENT_PREPROCESSING)
add_definitions_if(ENABLE_REDUCE_FOR_MANYCORE   -DARTED_REDUCE_FOR_MANYCORE)
add_definitions_if(ENABLE_LARGE_BLOCKING        -DARTED_LBLK)
add_definitions_if(USE_OLD_PROPAGATOR           -DARTED_USE_OLD_PROPAGATOR)

if (OPT_STENCIL)
  add_definitions(-DARTED_STENCIL_OPTIMIZED)

  add_definitions_if(STENCIL_WITH_C       -DARTED_STENCIL_WITH_C)
  add_definitions_if(DOMAIN_IS_POW2       -DARTED_DOMAIN_POWER_OF_TWO)
  add_definitions_if(ENABLE_EXPLICIT_VEC  -DARTED_EXPLICIT_VECTORIZATION)
  add_definitions_if(ENABLE_SWPREFETCH    -DARTED_ENABLE_SOFTWARE_PREFETCH)
  add_definitions_if(ENABLE_ARRAY_PADDING -DARTED_STENCIL_PADDING)
  add_definitions_if(ENABLE_LOOP_BLOCKING -DARTED_STENCIL_ENABLE_LOOP_BLOCKING)
else ()
  add_definitions(-DARTED_STENCIL_ORIGIN)
endif ()

add_definitions("${ADDITIONAL_MACRO}")

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/workaround.cmake)

# Numerical library
set(FILE_MATHLIB lapack)

### Build target
set(CORE_LIB src/core)
set(ARTED_LIB ARTED)
set(GCEED_LIB GCEED)

add_subdirectory(${ARTED_LIB})
add_subdirectory(${GCEED_LIB})
add_subdirectory(${CORE_LIB})

add_executable(${TARGET_NAME}
               main/main.f90
               main/read_input.f90
               main/setup_parallel.f90
               main/end_parallel.f90
               )

add_dependencies(${TARGET_NAME} ${SALMON_LINK_LIBRARIES})
target_link_libraries(${TARGET_NAME} ${SALMON_LINK_LIBRARIES} ${EXTERNAL_LIBS})

message(STATUS "Enabled macro")
get_directory_property(MACRO_LOG COMPILE_DEFINITIONS)
foreach(_macro ${MACRO_LOG})
  message(STATUS "  ${_macro}")
endforeach()