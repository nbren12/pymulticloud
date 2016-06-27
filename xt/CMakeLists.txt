cmake_minimum_required (VERSION 2.6)
project (HELLO Fortran)

set (CMAKE_Fortran_FLAGS_RELEASE "-O2")


SET ( SOURCES forcings.f90
    state_mod.f90
    central_scheme_correct.f90
    cmt_mod.f90
    multicloud_newform_relaxndstoch.f90
    stochastic.f90 util.f90
    duni.f
    )

add_subdirectory (ext)

add_executable ( mcnde ${SOURCES} )
target_link_libraries (mcnde minpack)
target_link_libraries (mcnde rk)
