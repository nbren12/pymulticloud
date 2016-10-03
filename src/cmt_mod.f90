module cmt_mod
  use util, only : searchsorted, cumsum
  implicit none
  public :: updatecmt, init_cmt, stochastic_cmt
  private

  logical :: stochastic_cmt = .false.

  ! storage arrays for adams-bashforth
  ! real(8), dimension(ntrunc, n) :: ko1 = 0d0, ko2 = 0d0
  real(8), dimension(:,:), allocatable :: ko1, ko2, tij

  real(8) :: pi, sqrt2
  real(8) :: taur, betaq, betau, qcref, hdref, duref, dumin, d0,dcmt,tauf

contains

  subroutine init_cmt(n, ntrunc)
    use param_mod
    integer, intent(in) :: n, ntrunc

    integer i, j, nz
    real(8) mult

    real(8) c, alpha_bar, L, T

    ! allocate arrays
    allocate(tij(0:16, ntrunc))
    allocate(ko1(ntrunc, n))
    allocate(ko2(ntrunc, n))
    
    pi = atan(1d0) * 4d0
    sqrt2 = sqrt(2.0d0)
    nz = size(tij,1)


    ! spectral to physical transformation matrix
    do j=lbound(tij,2),ubound(tij,2)
       do i=lbound(tij,1),ubound(tij,1)
          tij(i,j) = dcos(j * (i-1) * pi / (nz-1))
       end do
    end do


    !! initialize parameters


    zt=16.00*km !topospheric height
    zm = 5000.d0 !middle trop height
    zb=500.d0 !meters boundary layer depth
    zp = 8.d0 * km ! average height of penetrative clouds

    ! Scales
    c         = dsqrt(n2) *zt / pi
    alpha_bar = zm * n2 * theta0 / g ! temperature scale
    L         = dsqrt(c / beta)
    T         = L / c


    ! timescale of transitions
    taur = 8d0 * hour / t

    ! heating-based parameters
    betaq = 1d0/(10d0/alpha_bar /(day/t))
    hdref = 10d0/alpha_bar /(day/t)

    ! shear-based parameters
    mult = 1.0d0
    betau = 1d0/(10d0*mult/c)
    dumin = 5d0/c * mult
    duref = 20d0/c * mult

    ! damping and CMT strength
    dcmt = 1/(1.25d0 * day / t)
    d0   = 1d0 / (3d0 * day/t)



    ! namelist
    print *, 'd0=', d0, 'dcmt=', dcmt
    print * ,'betau, dumin, duref', betau, dumin, duref

  end subroutine init_cmt


  subroutine updatecmt(u, scmt, hd, hc, hs, dt)
    integer :: scmt(:)
    real(8) :: u(:,:), hd(:), hc(:), hs(:)
    real(8) :: dt
    intent(inout) u

    integer, parameter :: nz = 5

    ! work
    real(8) :: rate(2,2,size(u,2)), umid(size(u,2)), ulo(size(u,2))

    ! gillespie algorithm variables
    real(8) :: rsum(2), tau, t, runif
    integer :: action
    real(8) :: rands(size(u,2))
    real(8), dimension(size(u,1)) :: k1, k2 ,k3, k4
    integer i, j

    ! compute transition rates

    if (stochastic_cmt) then
       call trates(u, scmt, hd, rate)

       if (any(rate < 0d0) ) then
          print *, 'rate is less then zero... stop'
          stop -1
       end if

       call random_number(rands)

       ! Gillespie Algorithm
       do i=1,size(u,2)
          t = 0d0
          do
          call cumsum(rate(scmt(i), :, i), rsum)

          ! next jump time
          tau = - log(rands(i))/rsum(size(rsum))
          t = t +tau

          if (t < dt) then
             call random_number(runif)
             scmt(i) = searchsorted(rsum, runif * rsum(size(rsum)))
          else
             exit
          end if

          end do
       end do

       !! Cheaper algorithm
       ! do i=1,size(u,2)
       !    ! do switch if jump time is less then dt (this is a crude version of the Gillespie)
       !    if (rate(i) *dt < rands(i)) then
       !       if (trand < dt) then
       !          scmt(i) = 1 - scmt(i)
       !       end if
       !    end if
       ! end do

       !! Equilibrium
       ! do i=1,size(u,2)
       !    if ( rands(i) * (rate(0,1,i) + rate(1,0,i) ) < rate(1,0,i)) then
       !       scmt(i) = 0
       !    else
       !       scmt(i) = 1
       !    end if
       ! end do


    else
       scmt = 1
    end if

    do i=1,size(u,2)

       ! ! rk4
       ! call cmtforcing(u(:,i), hd(i), scmt(i), k1)
       ! call cmtforcing(u(:,i) + k1 * dt/2d0, hd(i), scmt(i), k2)
       ! call cmtforcing(u(:,i) + k2 * dt/2d0, hd(i), scmt(i), k3)
       ! call cmtforcing(u(:,i) + k3 * dt, hd(i), scmt(i), k4)
       ! u(:,i) = u(:,i) + dt/6d0 *(k1+ 2d0 * k2 + 2d0*k3 + k4)

       ! ! euler
       ! call cmtforcing(u(:,i), hd(i), scmt(i), k1)
       ! u(:,i) = u(:,i) + dt * k1

       ! ! trap
       ! call cmtforcing(u(:,i), hd(i), scmt(i), k1)
       ! call cmtforcing(u(:,i) + k1 *dt, hd(i), scmt(i), k2)
       ! u(:,i) = u(:,i) + dt/2d0 * (k1 + k2)

       ! ! AB-2
       ! call cmtforcing(u(:,i), hd(i), scmt(i), k1)
       ! call cmtforcing(uo1(:,i), hd(i), scmt(i), k2)
       ! u(:,i) = u(:,i) + dt/2d0 * (3d0 * k1 - k2)

       ! AB-3
       call cmtforcing(u(:,i), hd(i), scmt(i), k1)
       u(:,i) = u(:,i) + dt/12d0 * (23d0 * k1 - 16d0 * ko1(:,i) + 5d0*ko2(:,i))

       ko2(:,i) = ko1(:,i)
       ko1(:,i) = k1

    end do



  end subroutine updatecmt


  subroutine cmtforcing(u, hd, scmt, fu)
    real(8), intent(in) :: u(:)
    real(8), intent(in) :: hd
    integer, intent(in) :: scmt
    real(8), intent(out) :: fu(:)

    ! Work
    real(8) :: ulo, umid, kappa


    fu = 0d0

    if (scmt == 1) then
       fu = -u * d0
    else if (scmt == 2) then
       call dulowmid(matmul(tij,u), ulo, umid)
       if (umid * ulo < 0d0 ) then
          kappa = -(hd/hdref)**2 * umid * dcmt
       else
          kappa = 0d0
       end if

       fu(1) = kappa / sqrt2
       fu(3) = -kappa/sqrt2
    end if

  end subroutine cmtforcing

  subroutine trates(u, scmt, hd, rate)
    real(8), intent(in) :: u(:,:), hd(:)
    integer, intent(in) :: scmt(:)
    real(8), intent(out) :: rate(0:, 0:,:) ! probability rate of a switching

    ! Work
    real(8) :: uzg(lbound(tij,1):ubound(tij,1),size(u,2))
    real(8) :: rands(size(u,2))
    real(8) :: dumid, dulow
    integer :: i

    ! Compute grid-space transformation
    uzg = matmul(tij, u)
    rate = 0d0
    do i=1,size(uzg,2)
       call dulowmid(uzg(:,i), dulow, dumid)


       if (dulow > dumin) then
          rate(0,1,i) = exp(betau * abs(dulow)+ betaq * hd(i)) / taur
       else
          rate(0,1,i) = 0d0
       end if
       rate(1, 0, i) = exp(betau * (duref -abs(dulow)) + betaq * (hdref - hd(i)))/ taur
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
