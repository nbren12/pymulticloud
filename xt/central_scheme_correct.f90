
!     
!   NON OSCILLATORY CENTRAL SCHEME OF G-S JIANG AND EITAN TADMOR
!
SUBROUTINE Central_scheme(UC,DX,DT,N,Q_tld,alpha_tld,lmd_tld)
implicit none
INTEGER N, I,j
REAL*8 UC(5,-1:n+2), DX,DT,Q_tld,alpha_tld,lmd_tld

real*8 Ux(5,0:n+1),LMD, lmd_2,u_stag(5,-1:n+2) ,NLT
NLT=1.d0                  ! Nonlinearity set to zero.
LMD= DT/DX
LMD_2= LMD/2


!
!     get bound. cond. and slopes
!



call slopes(UC,UX,N)


!
!     second order central scheme
!
!
!     staggered cell averages reconstruction:   xi_stag(i,j):=xi(x_{i+1/2},y_{j+1/2})
!     
do J=1,5
    DO I=1,N
      u_stag(J,I)=(UC(J,I) + UC(J,I+1))/2 &
&           +(Ux(J,I)-Ux(J,I+1))/8
    ENDDO
ENDDO

!     
!     predictor: mid-time-step pointwise values at cell-center-grid-points (x_i)
!

do I=1,N


    UC(5,I) = UC(5,I) - LMD_2*( NLT*(UC(5,I)  + &
    &    Q_TLD)* UX(1,I) &
    &        + (alpha_tld*UC(5,I)*NLT + &
    &        lmd_tld*Q_TLD)* UX(2,I) &
    &      +  (UC(1,I)*NLT + alpha_tld*UC(2,I))*UX(5,I) &
    &        )

    UC(1,I) = UC(1,I) +  LMD_2*UX(3,I)
    UC(2,I) = UC(2,I) +  LMD_2*UX(4,I)
    UC(3,I) = UC(3,I) +  LMD_2*UX(1,I)
    UC(4,I) = UC(4,I) +  LMD_2*UX(2,I)/4
ENDDO





!
!     boundary conditions for mid-time-step cell-centered values
!


CALL PERIODIC_BC(UC, N)




!
!     Corrector: high resolution second order scheme
!

DO I=1,N

    u_stag(1,I)= u_stag(1,I) + LMD*(uC(3,I+1)-uC(3,I))
    u_stag(2,I)= u_stag(2,I) + LMD*(uC(4,I+1)-uC(4,I))
    u_stag(3,I)= u_stag(3,I) + LMD*(uC(1,I+1)-uC(1,I))
    u_stag(4,I)= u_stag(4,I) + LMD*(uC(2,I+1)-uC(2,I))/4

    u_stag(5,I)= u_stag(5,I) - LMD*( &
    &        UC(5,I+1)*(UC(1,I+1) + ALPHA_TLD*UC(2,I+1))*NLT + &
    &        Q_TLD*(UC(1,I+1) + LMD_TLD*UC(2,I+1)) &
    &        -UC(5,I)*(UC(1,I)+ ALPHA_TLD*UC(2,I))*NLT &
    &        - Q_TLD*(UC(1,I) + LMD_TLD*UC(2,I)))

ENDDO

!
!     Boundray conditions and slopes for staggered values
!


call slopes(U_STAG,UX,N)

!
!     Centered cell averaging:
!

do J=1,5
    DO I=1,N
      uC(J,I)=(U_STAG(J,I-1) + U_STAG(J,I))/2 &
      &           +(Ux(J,I-1)-Ux(J,I))/8
    ENDDO
ENDDO

RETURN
END

SUBROUTINE PERIODIC_BC(UC, N)
implicit none
INTEGER J,N
REAL*8 UC(5,-1:n+2)
DO J=1,5
    UC(J,-1) = UC(J,N-1)
    UC(J,0) =  UC(J,N)
    UC(J,N+1) = UC(J,1)
    UC(J,N+2) = UC(J,2)
ENDDO
RETURN
END


subroutine slopes(uc,ux,N)
implicit none
integer N,I,J
!      parameter(nmax=256,mmax=200)
real*8 uc(5,-1:n+2)
real*8 ux(5,0:n+1),tht,bminmod
parameter(tht=1.5d0)      !tht must be between 1 to 2.

logical limiter

parameter(limiter=.true.)

!
!     periodic boundary conditions in x
!

call  PERIODIC_BC(UC, N)



!
!    High resolution slopes calculation
!

do j=1,5
    do i=0,N+1
      if(limiter)then
          ux(j,i)= bminmod(tht*(uc(j,i+1)-uc(j,i)), &
          &               (uc(j,i+1)-uc(j,i-1))/2,tht*(uc(j,i)-uc(j,i-1)))


      else
          ux(j,i)= (uc(j,i+1)-uc(j,i-1))/2


      endif
    enddo
enddo




return
end



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
end



