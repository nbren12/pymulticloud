subroutine cmt_wrapper(u, scmt, hd, hc, hs, n, ntrunc, dt) bind(c)
  use, intrinsic :: iso_c_binding
  use cmt_mod
  implicit none

  integer(c_int), intent(in)    :: n, ntrunc
  real(c_double), intent(inout) :: u(ntrunc, n)
  real(c_double), intent(in)    :: hd(n), hc(n), hs(n)
  integer(c_int), intent(inout)    :: scmt(n)
  real(c_double), intent(in)    :: dt

  logical :: first_call = .true.

  if (first_call) then
     call init_cmt(n, ntrunc)
  end if

  call updatecmt(u, scmt, hd, hc, hs, dt)
end subroutine cmt_wrapper
