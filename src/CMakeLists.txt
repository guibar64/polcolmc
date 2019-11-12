
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")


set(CORE_SRC
  core/random2.f90
  core/readparam.f90
  core/geom.f90
  core/celldec.f90
  core/inputp.f90
  core/state.f90
  core/misc.f90
)

set(NRG_SRC
  nrg/vscoul2pc.F90
  nrg/gauss.f90
  nrg/none.F90
  nrg/yukawa.F90
  nrg/yukawa2.F90
  nrg/fennel.F90
  nrg/nrg.F90
)

set(SIMU_SRC
  moves.f90
  simulation.f90
  spemoves.f90
  analysis/sinft.F90
  analysis/correls.F90
  analysis/struct_fact.f90
  analysis/spline.f90
  analysis/agr.f90
  analysis/conftemp.f90
  analysis/fluct_dens.f90
  analysis/sdeq_comp.f90
  analysis/spherharm.f90
  analysis/bondorder.f90
  analysis/zdens.f90
  analysis/chempot.f90
  analysis/swapmap.f90
  analysis/statcell.f90
  analysis/ergo.f90
  analysis/histograms.f90
  tools/trr.F90
  tools/tools.f90
  sgc.f90
)

add_library(pcmc ${SIMU_SRC} ${NRG_SRC} ${CORE_SRC})
set_property(TARGET pcmc APPEND_STRING PROPERTY Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
target_compile_definitions(pcmc PUBLIC FOURIER_TRANSFORM=${FOURIER_TRANSFORM})

include(TestBigEndian)
TEST_BIG_ENDIAN(BIG_ENDIAN)
if(BIG_ENDIAN)
target_compile_definitions(pcmc PUBLIC "POLCOLMC_BIG_ENDIAN=1")
endif()

# Check support for CONVERT specifier in open() stmt (for trr files)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
message(STATUS "CONVERT specifier is supported")
set(CONVERSPEC_SUPPORT TRUE)
else()
endif()

if(TRR_SUPPORT STREQUAL "FALSE" OR NOT CONVERSPEC_SUPPORT)
target_compile_definitions(pcmc PUBLIC "NO_CONVERSPEC_SUPPORT")
message(STATUS "Will build without support for trr files.")
else()
message(STATUS "Will build with support for trr files.")
endif()


if(OpenMP_FOUND)
target_compile_options(pcmc PUBLIC "${OpenMP_Fortran_FLAGS}")
endif()

