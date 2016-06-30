module cmt_mod
  use state_mod, only: ntrunc
  implicit none
  public :: updatecmt, init_cmt
  private

  real(8) :: pi
  real(8) :: taur, betaq, betau, qcref, hdref, duref, dumin, d1, d2, tf

  real(8) :: tij(0:16, ntrunc)
contains

  subroutine init_cmt
    use param, only: t, hour, day, c, alpha_bar
    integer i, j, nz


    pi = atan(1d0) * 4d0
    nz = size(tij,1)

    ! spectral to physical transformation matrix
    do j=lbound(tij,2),ubound(tij,2)
       do i=lbound(tij,1),ubound(tij,1)
          tij(i,j) = dcos(j * (i-1) * pi / (nz-1))
       end do
    end do


    !! initialize parameters

    ! timescale of transitions
    taur = 8d0 * hour / t

    ! heating-based parameters
    betaq = 1d0/(10d0/alpha_bar /(day/t))
    hdref = 10d0/alpha_bar /(day/t)

    ! shear-based parameters
    betau = 1d0/(10d0/c)
    dumin = 5d0/c
    duref = 20d0/c


  end subroutine init_cmt


  subroutine updatecmt(u, scmt, hd, hc, hs, dt)
    integer :: scmt(:)
    real(8) :: u(:,:), hd(:), hc(:), hs(:)
    real(8) :: dt

    integer, parameter :: nz = 5

    ! work
    real(8) :: rate(size(u,2))
    real(8) :: rands(size(u,2))
    real(8) :: pr
    integer i, j, trand


    ! compute transition rates
    call trates(u, scmt, hd, rate)

    if (any(rate < 0d0) ) then
       print *, 'rate is less then zero... stop'
       stop -1
    end if

    call random_number(rands)

    do i=1,size(u,2)

       ! do switch if jump time is less then dt (this is a crude version of the Gillespie)
       if (rate(i) > 1d-10) then
          trand = -log(rands(i))/rate(i)
          print *, rands(i)
          if (trand < dt) then
             scmt(i) = 1 - scmt(i)
          end if
       end if
    end do

  end subroutine updatecmt

  subroutine trates(u, scmt, hd, rate)
    real(8), intent(in) :: u(:,:), hd(:)
    integer, intent(in) :: scmt(:)
    real(8), intent(out) :: rate(:) ! probability rate of a switching

    ! Work
    real(8) :: uzg(lbound(tij,1):ubound(tij,1),size(u,2))
    real(8) :: dumid, dulow
    integer :: i

    ! Compute grid-space transformation
    uzg = matmul(tij, u)

    do i=1,size(uzg,2)
       call dulowmid(uzg(:,i), dulow, dumid)

       if (scmt(i) == 0) then
          if (dulow > dumin) then
             rate(i) = exp(betau * abs(dulow)+ betaq * hd(i)) / taur
          else
             rate(i) = 0d0
          end if
       else if (scmt(i) == 1) then
          rate(i) = exp(betau * (duref -abs(dulow)) + betaq * (hdref - hd(i)))/ taur
       else
          print *, 'Invalid value: scmt(', i, ') =', scmt(i)
          stop -1
       end if

    end do

  end subroutine trates

  subroutine dulowmid(uz, dulow, dumid)
    real(8), intent(in) :: uz(0:)
    real(8), intent(out) :: dulow, dumid

    integer i, ilow, ihi, lowstart, histart

    lowstart = 1
    histart  = 7

    ilow = maxloc(abs(uz(lowstart:7) - uz(0)), dim=1) + lowstart -1
    dulow = uz(ilow) - uz(0)

    ihi = maxloc(abs(uz(histart:13) - uz(ilow)), dim=1) + histart -1
    dumid = uz(ihi) - uz(ilow)

  end subroutine dulowmid

end module cmt_mod
