cmake_minimum_required (VERSION 2.6)
project (HELLO Fortran)

set (CMAKE_Fortran_FLAGS_RELEASE "-O2")

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS_DEBUG "-check bounds -traceback -g")
endif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")



SET ( SOURCES
  param_mod.f90
  state_mod.f90
  cmt_mod.f90
  forcings.f90
  stochastic.f90 central_scheme_correct.f90
  util.f90
  duni.f
  multicloud_newform_relaxndstoch.f90
  )

add_subdirectory (ext)

add_executable ( mcnde ${SOURCES} )
target_link_libraries (mcnde minpack)
target_link_libraries (mcnde rk)
