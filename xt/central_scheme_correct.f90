!
!   non oscillatory central scheme of g-s jiang and eitan tadmor
!
module nonlinear_module
  implicit none

  logical, parameter :: toggle_nonlinear = .false.
  integer, parameter :: ntrunc = 2
contains

  subroutine vertical_advection_temp(w,t,f)

    real(8), intent(in) :: w(:), t(:)
    real(8), intent(out) :: f(:)

    integer i,m,L
    L = size(f)

    do m=1,L
       f(m) = 0

       do i=1,m-1
          f(m) = f(m) + t(i) * w(m-i)
       end do

       do i=1,L-m
          f(m) = f(m) - w(i) * t(i+m) - t(i) * w(i+m)
       end do

       f(m) = f(m) * m / dsqrt(2.0d0)
    end do

  end subroutine vertical_advection_temp

  subroutine vertical_advection_u(w,u,f)

    real(8), intent(in) :: w(:), u(:)
    real(8), intent(out) :: f(:)

    integer i,m,L
    L = size(f)

    do m=1,L
       f(m) = 0

       do i=1,m-1
          f(m) = f(m) + u(i) * w(m-i)
       end do

       do i=1,L-m
          f(m) = f(m) - w(i) * u(i+m) + u(i) * w(i+m)
       end do

       f(m) = f(m) * m / dsqrt(2.0d0)
    end do

  end subroutine vertical_advection_u

  subroutine vertical_advection_tendency(uc, f, dx,dt)
    real(8) uc(:,-1:)
    real(8) dx, dt
    integer i, j, n

    real(8) f(:,:), w(ntrunc), temp(ntrunc)

    n  = ubound(uc, 2) -2

    call periodic_bc(uc, 2)

    if (toggle_nonlinear) then
       do i=1,n

          ! read in from state array
          do j=1,ntrunc
             w(j) = (uc(j, i-1) - uc(j, i+1))/dx
             temp(j) = uc(ntrunc+j,i)
          end do

          call vertical_advection_u(w, uc(1:ntrunc,i), f(1:ntrunc, i))
          call vertical_advection_temp(w, temp, f(ntrunc+1:2*ntrunc, i))
       end do
    end if

  end subroutine vertical_advection_tendency

  subroutine vertical_advection_driver(uc, dx, dt, n)
    real(8) uc(5,-1:n+2), dx,dt
    integer n

    real(8) f(size(uc,1), size(uc,2), 2)

    if (toggle_nonlinear) then
       call vertical_advection_tendency(uc, f(:,:,1), dx, dt)

       uc = uc + dt * f(:,:,1)

       call vertical_advection_tendency(uc, f(:,:,2), dx, dt)
       uc = uc + dt/2d0 * (f(:,:,2) - f(:,:,1))
    end if

  end subroutine vertical_advection_driver

  subroutine central_scheme(uc,dx,dt,n,q_tld,alpha_tld,lmd_tld)
    implicit none
    integer n, i,j
    real*8 uc(5,-1:n+2), dx,dt,q_tld,alpha_tld,lmd_tld, fc(5,-1:n+2)

    real*8 ux(5,0:n+1),lmd, lmd_2,u_stag(5,-1:n+2) ,nlt
    nlt=1.d0                  ! nonlinearity set to zero.
    lmd= dt/dx
    lmd_2= lmd/2


    !
    !     get bound. cond. and slopes
    !



    call slopes(uc,ux)


    !
    !     second order central scheme
    !
    !
    !     staggered cell averages reconstruction:   xi_stag(i,j):=xi(x_{i+1/2},y_{j+1/2})
    !
    do j=1,5
       do i=1,n
          u_stag(j,i)=(uc(j,i) + uc(j,i+1))/2 &
               &           +(ux(j,i)-ux(j,i+1))/8
       enddo
    enddo

    !
    !     predictor: mid-time-step pointwise values at cell-center-grid-points (x_i)
    !


    do i=-1,n+2
       call flux(uc(:,i), fc(:,i), q_tld,alpha_tld,lmd_tld)
    end do
    call slopes(fc, ux)

    do i=1,n
       do j =1,5
          uc(j,i) = uc(j,i) - lmd_2 * ux(j,i)
       end do
    enddo





    !
    !     boundary conditions for mid-time-step cell-centered values
    !


    call periodic_bc(uc, 2)

    do i=-1,n+2
       call flux(uc(:,i), fc(:,i), q_tld,alpha_tld,lmd_tld)
    end do




    !
    !     corrector: high resolution second order scheme
    !

    do i=1,n
       do j=1,5

       u_stag(j,i)= u_stag(j,i) - lmd*(fc(j,i+1)-fc(j,i))
       end do
    enddo

    !
    !     boundray conditions and slopes for staggered values
    !


    call slopes(u_stag,ux)

    !
    !     centered cell averaging:
    !

    do j=1,5
       do i=1,n
          uc(j,i)=(u_stag(j,i-1) + u_stag(j,i))/2 &
               &           +(ux(j,i-1)-ux(j,i))/8
       enddo
    enddo

    return
  end subroutine central_scheme

  subroutine flux(u, f, q_tld,alpha_tld,lmd_tld)

    real(8), intent(in) :: u(:)
    real(8), intent(out):: f(:)
    real*8 q_tld,alpha_tld,lmd_tld


    integer m, i, L
    real(8) theta(ntrunc), q
    real(8) ftheta(ntrunc), fq


    L = ntrunc


    do m=1,ntrunc
       theta(m) = u(m+L)
    end do

    q = u(2*L+1)

    f = 0.0d0
    ftheta = 0.0d0
    fq = 0.0d0

    do m=1,ntrunc

       if (toggle_nonlinear) then

          ! velocity
          do i=1,m-1
             f(m) = f(m) + u(i) * u(m-i)
          end do

          do i=1,L-m
             f(m) = f(m) + 2 * u(i) * u(i+m)
          end do

          f(m) = f(m)/ dsqrt(2.0d0)

          ! temperature
          do i=1,m-1
             ftheta(m) = ftheta(m) + u(i) * theta(m-i)
          end do

          do i=1,L-m
             ftheta(m) = ftheta(m) +  u(i) * theta(i+m) - u(i+m) * theta(i)
          end do

          ftheta(m) = ftheta(m)/ dsqrt(2.0d0)

       end if


       ! Linear terms

       ! Pressure gradient
       f(m) = f(m) - theta(m) / m

       ! convergence
       ftheta(m) = ftheta(m) - u(m) / m


    end do

    ! moisture equation
    fq  = q_tld * (u(1) + lmd_tld * u(2)) + q *(u(1) + alpha_tld * u(2))


    do m=1,ntrunc
       f(L+m) = ftheta(m)
    end do

    f(2*L+1) = fq

  end subroutine flux

  subroutine periodic_bc(phi, g)
    real(8) :: phi(:,:)
    integer n,g,i,j
    n = size(phi,2)


    do i=1,g
       do j=1,size(phi,1)
          phi(j,i) = phi(j,n-g-g+1)
          phi(j,n-g+i) = phi(j,g+i)
       end do
    end do

  end subroutine periodic_bc

  subroutine slopes(uc,ux)
    implicit none
    integer n,i,j
    !      parameter(nmax=256,mmax=200)
    real*8 uc(:,-1:)
    real*8 ux(:,0:),tht
    parameter(tht=1.5d0)      !tht must be between 1 to 2.

    logical limiter

    parameter(limiter=.true.)

    !
    !     periodic boundary conditions in x
    !

    call  periodic_bc(uc, 2)



    !
    !    high resolution slopes calculation
    !

    do j=1,size(uc,1)
       do i=0,ubound(ux,2)
          if(limiter)then
             ux(j,i)= bminmod(tht*(uc(j,i+1)-uc(j,i)), &
                  &               (uc(j,i+1)-uc(j,i-1))/2,tht*(uc(j,i)-uc(j,i-1)))


          else
             ux(j,i)= (uc(j,i+1)-uc(j,i-1))/2


          endif
       enddo
    enddo




    return
  end subroutine slopes



  real*8 function bminmod(a,b,c)
    real*8 a,b,c
    if(a.gt.0.d0.and.b.gt.0.d0.and.c.gt.0.d0) then

       bminmod= dmin1(a,b,c)
    elseif(a.lt.0.d0.and.b.lt.0.d0.and.c.lt.0.d0) then
       bminmod= dmax1(a,b,c)
    else
       bminmod=0.d0
    endif
    return
  end function bminmod



end module nonlinear_module
