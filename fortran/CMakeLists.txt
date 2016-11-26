cmake_minimum_required (VERSION 2.6)
project (HELLO Fortran)

set (CMAKE_Fortran_FLAGS_RELEASE "-O2")

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS_DEBUG "-check bounds -traceback -g")
endif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")



SET ( SOURCES
  multicloud_mod.f90
  param_mod.f90
  state_mod.f90
  cmt_mod.f90
  central_scheme_correct.f90
  util.f90
  duni.f)

add_subdirectory (ext)
add_executable ( mcnde ${SOURCES} multicloud_newform_relaxndstoch.f90)

add_executable ( column ${SOURCES} column_wrapper.f90)
add_library(multicloud SHARED multicloud_wrapper.f90 cmt_wrapper.f90 cmt_mod.f90 multicloud_mod.f90 param_mod.f90 util.f90 duni.f )

target_link_libraries (mcnde minpack)
target_link_libraries (multicloud minpackdl)
target_link_libraries (column minpack)
