subroutine multicloud_wrapper(fc, fd, fs, u1, u2, t1, t2, teb, q, hs, n, dt,&
     dx, time, tebst) bind(c)
  use, intrinsic :: iso_c_binding
  use util
  use param_mod
  use multicloud_mod, only: init_multicloud, updatehcds, range_kuttas, get_eqcloudfrac, T, L, alpha_bar, c, ud, theta_ebs_m_theta_eb

  IMPLICIT NONE
  INTEGER :: i, j, niter, iter,k,nout, iout,inrg
  logical :: first_call = .true.


  !     PROGNOSTIC VARIABLES
  integer(c_int), intent(in) :: n
  real(c_double), intent(in) :: tebst(n)
  real(c_double), intent(inout) ::  u1(n),u2(n),t1(n),t2(n),q(n),teb(n),hs(n),&
       fc(n), fd(n), fs(n)

  real(c_double), intent(in) :: dt, dx, time

  real(8) :: hds(n), hd(n), hc(n)


  real(8) fceq, fdeq, fseq

  if (first_call) then
     call init_param()
     call init_multicloud
     ! call get_eqcloudfrac(fceq, fdeq, fseq)
  end if

  CALL updatehcds(fc,fd,fs,u1, u2, t1,t2,teb,q,hds,hc,hd  &
       ,n, dx, 2.d0*dt*t/(hour) )
  CALL range_kuttas(u1,u2,t1,t2,teb,q,hs,hc,hd,n,&
       dt,tebst,time,hds)
  first_call = .false.
end subroutine multicloud_wrapper

subroutine multicloud_eq(fceq, fdeq, fseq) bind(c)
  use, intrinsic :: iso_c_binding
  use util
  use param_mod
  use multicloud_mod, only: init_multicloud, updatehcds, range_kuttas, get_eqcloudfrac, T, L, alpha_bar, c, ud, theta_ebs_m_theta_eb

  IMPLICIT NONE
  INTEGER :: i, j, niter, iter,k,nout, iout,inrg
  logical :: first_call = .true.


  !     PROGNOSTIC VARIABLES
  real(c_double), intent(out) ::  fceq, fdeq, fseq


  call init_param()
  call init_multicloud
  call get_eqcloudfrac(fceq, fdeq, fseq)

end subroutine multicloud_eq
