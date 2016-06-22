!     
!   non oscillatory central scheme of g-s jiang and eitan tadmor
!
module nonlinear_module
  implicit none
contains
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


    call flux(uc, fc,q_tld,alpha_tld,lmd_tld)
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
    call flux(uc, fc,q_tld,alpha_tld,lmd_tld)




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

  subroutine flux(uc,fc,q_tld,alpha_tld,lmd_tld)
    real(8), intent(in) :: uc(:,:)
    real(8), intent(out) :: fc(:,:)
    real*8 q_tld,alpha_tld,lmd_tld
    integer i

    do i =1,size(uc,2)
       fc(1,i) = -uc(3,i)
       fc(2,i) = -uc(4,i)
       fc(3,i) = -uc(1,i)
       fc(4,i) = -uc(2,i)/4
       fc(5,i) = q_tld * (uc(1,i) + lmd_tld * uc(2, i)) + uc(5,i) * (uc(1,i) + alpha_tld * uc(2,i))
    end do
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
