cmake_minimum_required (VERSION 2.8)
project (PolyMCsim Fortran C)

set(PROJ_CONTACT "guillaume.bareigts@ntymail.com")

set(PCMSIM_SUMMARY "PolyMCsim is a PhD project to perform Monte-Carlo simulations of polydisperse colloidal dispersions using effective pair potentials (Yukawa).")

set(POLCOLMC_VERSION_MAJOR 0)
set(POLCOLMC_VERSION_MINOR 7)
set(POLCOLMC_VERSION_PATCH 0)

cmake_host_system_information(RESULT HOST_NAME QUERY HOSTNAME)

# git info
execute_process(                                                  
  COMMAND "git" "rev-parse" "HEAD"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_SHA1
  ERROR_QUIET                                                           
  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(
  COMMAND "git" "rev-parse" "--abbrev-ref" "HEAD"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  ERROR_QUIET                                                           
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if(GIT_SHA1 STREQUAL "")                       
  set(POLCOLMC_VERSION_LONG "${POLCOLMC_VERSION_MAJOR}.${POLCOLMC_VERSION_MINOR}.${POLCOLMC_VERSION_PATCH}")
else()
  set(POLCOLMC_VERSION_LONG "${POLCOLMC_VERSION_MAJOR}.${POLCOLMC_VERSION_MINOR}.${POLCOLMC_VERSION_PATCH}-git${GIT_SHA1}")
endif() 

set(POLCOLMC_VERSION "${POLCOLMC_VERSION_MAJOR}.${POLCOLMC_VERSION_MINOR}.${POLCOLMC_VERSION_PATCH}")

# Versionned prefix
if(VERSION_PREFIX STREQUAL "LONG")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}-${POLCOLMC_VERSION_LONG}")
elseif(VERSION_PREFIX STREQUAL "SHORT")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}-${POLCOLMC_VERSION}")
else()
  set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
endif()
message(STATUS "Prefix: ${CMAKE_INSTALL_PREFIX}")

if(FOURIER_TRANSFORM STREQUAL "FFTW3")
  # try to found fftw3
  find_path(FFTW_INCLUDE NAMES "fftw3.h")
  if(FFTW_INCLUDE STREQUAL "FFTW_INCLUDE-NOTFOUND")
    message(FATAL_ERROR "fftw3.h : not found")
  endif()
  find_library(FFTW_LIB NAMES fftw3)
  if(FFTW_LIB STREQUAL "FFTW_LIB-NOTFOUND")
    message(FATAL_ERROR "fftw3 library : not found")
  else()
    message(STATUS "Will use ffwt3 lib: ${FFTW_LIB}")
  endif()
else()
  set(FOURIER_TRANSFORM "SLOW")
  set(FFTW_INCLUDE "")
  set(FFT_LIBS "")
endif()
set(FFT_INCLUDE "${FFTW_INCLUDE}")
set(FFT_LIBS "${FFTW_LIB}")
include_directories(SYSTEM ${FFT_INCLUDE})

if(OPENMP)
find_package(OpenMP)
if(OpenMP_FOUND)
message(STATUS "Using OpenMP")
endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

# pcmc lib
add_subdirectory(src)

if(NOT POLCOLMC_EXE)
  set(POLCOLMC_EXE "polcolmc")
endif(NOT POLCOLMC_EXE)

add_executable(${POLCOLMC_EXE} src/main.F90)
set_property(TARGET ${POLCOLMC_EXE} APPEND_STRING PROPERTY Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
target_compile_definitions(${POLCOLMC_EXE} PUBLIC __COMMIT_HASH__="${GIT_SHA1}" FOURIER_TRANSFORM=${FOURIER_TRANSFORM})

target_link_libraries(${POLCOLMC_EXE} pcmc ${FFT_LIBS})
if(OpenMP_FOUND)
target_compile_options(${POLCOLMC_EXE} PUBLIC "${OpenMP_Fortran_FLAGS}")
target_link_libraries(${POLCOLMC_EXE} "${OpenMP_Fortran_LIBRARIES}")
endif()

install(TARGETS ${POLCOLMC_EXE} DESTINATION bin)

# tests
enable_testing()
add_subdirectory(tests)

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
  "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_VERSION_MAJOR "${POLCOLMC_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${POLCOLMC_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${POLCOLMC_VERSION_PATCH}")
set(CPACK_PACKAGE_CONTACT "${PROJ_CONTACT}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${POLCOLMC_SUMMARY}")
set(CPACK_PACKAGE_NAME "polcolmc")
include (CPack)
