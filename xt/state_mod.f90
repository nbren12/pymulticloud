module state_mod
  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: n = 1040



  ! stochatic arrays
  integer scmt(n) ! CMT state

  ! cloud fraction
  real(dp) :: fcls(n), fdls(n), fsls(n)


contains

  subroutine initialize_scmt()
    scmt = 0
  end subroutine initialize_scmt
end module state_mod
