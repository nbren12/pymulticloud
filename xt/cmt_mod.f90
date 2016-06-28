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

    integer, parameter :: nz = 5
    logical :: first_call = .true.

    ! work
    real(8) :: uzg(nz,size(u,2)), duh(size(u,2)), dul(size(u,2))
    integer :: zs(size(u,2))
    real(8) :: R(1:3, 1:3)
    real(8), save :: zij(nz,nz)
    real(8) :: pi
    integer i, j

    if (first_call) then

       pi = datan(1.0d0) * 4.d0
       print *, 'pi = ', pi


       do i=1,nz
          do j=1,nz
             zij(i,j) = dcos(j * (i-1) * pi / (nz-1))
          end do
       end do

       first_call = .false.
    end if

    ! Compute grid-space transformation
    uzg = matmul(zij, u)


  end subroutine updatecmt

end module cmt_mod
