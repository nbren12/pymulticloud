module cmt_mod
  implicit none
  public :: updatecmt
  private

  real(8) :: taur, betalam, betaq, betau, qcref, qdref, duref, dumin, d1, d2, tf
contains
  subroutine updatecmt(u, scmt, hd, hc, hs, dt)
    integer :: scmt(:) 
    real(8) :: u(:,:), hd(:), hc(:), hs(:)
    real(8) :: dt

    integer, parameter :: nz = 10

    ! work
    real(8) :: uzg(nz,size(u,2)), duh(size(u,2)), dul(size(u,2))
    real(8) :: R(1:3, 1:3)
    integer :: zs(size(u,2))
    integer i


    do i=1,nz
       uzg(i, )
    end do

  end subroutine updatecmt

end module cmt_mod
