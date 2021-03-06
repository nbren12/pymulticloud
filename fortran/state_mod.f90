module state_mod
  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: n = 200, ntrunc=2



  ! stochatic arrays
  integer scmt(n) ! CMT state

  ! cloud fraction
  real(dp) :: fcls(n), fdls(n), fsls(n)

  ! lmd
  real(dp) :: lmd(n)
end module state_mod
