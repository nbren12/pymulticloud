module state_mod
  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: n = 520, ntrunc=4



  ! stochatic arrays
  integer scmt(n) ! CMT state

  ! cloud fraction
  real(dp) :: fcls(n), fdls(n), fsls(n)
end module state_mod
