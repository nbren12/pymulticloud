module param
real*8 nstochgl

real*8 a,b,a0,a1,a2,a0p,a1p,a2p,alpha2,alpha3,  &
    lambdas,alpha_s,alpha_c,xis,xic,mu,ud,thetap,thetam,zb,zt
real*8 tau_conv,tau_e, tau_s,tau_c,tau_r, tau_d, t, sdc,dtype

REAL*8 tau01,tau02,tau10,tau12,tau20,tau30,tau23, r23value, times
REAL*8 capebar, cape0,dry,pbar, qbar,fceq,fdeq,fseq, rstoch

real*8 theta_eb_elbar,deltac,deltac1

real*8 qr01,qr02,lambdabar,qcbar,dbar,  &
    theta_eb_m_theta_em,m0,theta_ebs_m_theta_eb,moist0,alpha4

COMMON /rce_values/qr01,qr02,lambdabar,qcbar,dbar,  &
    theta_eb_m_theta_em,m0,theta_ebs_m_theta_eb,moist0,alpha4

COMMON/const_par/a,b,a0,a1,a2,a0p,a1p,a2p,alpha2,alpha3,  &
    lambdas,alpha_s,alpha_c,xis,xic,mu,ud,thetap,thetam,zb,zt

COMMON /char_times/tau_conv,tau_e, tau_s,tau_c,tau_r,t, sdc,dtype

COMMON /tauscales/ capebar, cape0,dry,pbar, qbar,fceq,fdeq,fseq
COMMON /tauscales/  tau01,tau02,tau10,tau12, rstoch
COMMON /tauscales/  tau20,tau30,tau23, r23value, times,nstochgl
COMMON /congd/ theta_eb_elbar,deltac,deltac1

contains

subroutine setup_param()
! Calculate values that depend on others
    a=(1.D0-Lambdas)/(thetap-thetam)
    b=Lambdas-a*thetam            !linear fit coefficients of LAMBDA

end subroutine
end module param

module forcings
use param
implicit none
!use const_par
!use rce_values



contains

subroutine DU(u, F, T, n)
! Order of u is u_1, u_2, t_1, t_2, teb, q, hs, hc,thteb_s
    integer n
    real *8 U(n), F(n)
    real *8 T
    intent(out) :: F

    real*8 q, theta1, theta2, theta_eb, hs, hc,thteb_st
    real *8 pr0, pr1, tht_eb_tht_em, lambda
    real*8 Hd, hsbar, hcbar, D, Pr
    real *8 pi, two_sqrt2

    theta1 =u(3)
    theta2 =   u(4)
    theta_eb = u(5)
    q        = u(6)
    hs =u(7)
    hc =u(8)
    thteb_st = u(9)

    TWO_SQRT2 = 2.0d0 * dsqrt(2.0d0)
    pi = 4*datan(1.d0)

    !--------Velocity Drag---------------------------------------



    !----------Deep Heating--------------------------------------

    ! Calculate Q_d
     PR0 = QcBar + (  A1*THETA_EB + A2*Q - &
         A0*(THETA1+ ALPHA2*THETA2))/TAU_CONV

     pr0 = dmax1(pr0,0.d0)

    ! Calculate Q_c
     PR1 = QcBar + (  A1p*THETA_EB + A2p*Q - &
         A0p*(THETA1+ ALPHA2*THETA2))/TAU_CONV
     PR1=dmax1(PR1,0.d0)

     ! this adds back the RCE quantity
     THT_EB_THT_EM = THETA_EB_M_THETA_EM + THETA_EB  &
         - (TWO_SQRT2/pi)*(THETA1 + ALPHA3*THETA2)- Q

     lambda = moisture_switch(THT_EB_THT_EM)

     ! Deep heating
     Hd =(1.D0 - LAMBDA)/(1.D0-LAMBDAS)*PR0

     !--------Downdrafts-----------------------------------------

     hsbar = alpha_s*(1.d0 -Lambdabar)*QcBar/(1.d0-Lambdas)
     hcbar = alpha_c*(Lambdabar-Lambdas)*QcBar/(1.d0-Lambdas)
     D = DMAX1(0.D0,(1.D0 + MU*( hsbar -hcbar +HS- HC)) )*m0 &
        *THT_EB_THT_EM


     !--------Precip---------------------------------------------

     PR=  Hd + xic*(hcbar+Hc) + xis*(hsbar + hs )


     !--------Calculate Source Quantities------------------------

     F(1) = -UD * u(1)
     F(2) = -UD * u(2)
     F(3) = PR - QR01 - THETA1/TAU_R
     ! This was an interesting detail in Boualem's code (no qr02, hsbar or hcbar
     ! tricky thing where boualem used the fact taht hcbar + hsbar - QR02 is
     ! zero. Not a bug per se, but very confusing when reading the code.
     F(4) = ( HC + hcbar ) - (HS + hsbar )- QR02 -theta2 / tau_r
     F(5) = (thteb_st- THETA_EB)/TAU_E - (D -DBAR)/ZB
     F(6) = D/ZT - TWO_SQRT2*PR/pi
     F(7) = (ALPHA_S* Hd - hsbar -HS)/TAU_S
     F(8) = (ALPHA_C*(LAMBDA-LAMBDAS)/(1.d0-LAMBDAS)*PR1&
             - hcbar - HC)/TAU_C
     F(9) = 0.0d0


end subroutine

subroutine step_source(u1,u2,theta1,theta2,theta_eb,q,hs,hc,hd, &
    n,DT, thteb_st)
    integer i,n
    real*8 u1(N),u2(N),theta1(N),theta2(N),hd(N), &
      theta_eb(N),q(N),hs(N),hc(N),thteb_st(n), DT
    real*8 u(9)

    do i = 1,n
        u = (/ u1(i), u2(i), theta1(i), theta2(i), theta_eb(i), q(i),&
        hs(i), hc(i), thteb_st(i)/)

!        call rk4(DU, u, 0.0d0, dt, 9, 1)
    end do


end subroutine


double precision function moisture_switch(x)
    real *8 x, lambda

    IF(x.GT.THETAP)THEN
       LAMBDA=1.D0
        !print*,'dry RCE. stoped'
        !stop
    ELSEIF (x.GT.THETAM)THEN
       LAMBDA = A* x + B
        !print*,'mixed RCE. stoped'
        !stop
    ELSE
       LAMBDA = LAMBDAS
    ENDIF
    moisture_switch = LAMBDA
end function moisture_switch

SUBROUTINE range_kuttas(u1,u2,theta1,theta2,theta_eb,q,hs,hc,hd,  &
        n,dt, thteb_st,timel,hds,s1)
use param
IMPLICIT NONE


REAL*8, INTENT(IN OUT)                   :: u1(n)
REAL*8, INTENT(IN OUT)                   :: u2(n)
REAL*8, INTENT(IN OUT)                   :: theta1(n)
REAL*8, INTENT(IN OUT)                   :: theta2(n)
REAL*8, INTENT(IN OUT)                   :: theta_eb(n)
REAL*8, INTENT(IN OUT)                   :: q(n)
REAL*8, INTENT(IN OUT)                   :: hs(n)
REAL*8, INTENT(IN)                       :: hc(n)
REAL*8, INTENT(IN)                       :: hd(n)
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN)                       :: dt
REAL*8, INTENT(IN)                       :: thteb_st(n)
REAL*8, INTENT(IN)                       :: timel
REAL*8, INTENT(IN OUT)                   :: hds(n)
REAL*8, INTENT(IN)                       :: s1
INTEGER :: i

REAL*8 ftht1,ftht2,fthteb,fq,fhs,fhc,two_sqrt2,xgtemp
REAL*8 pr0,pr1,pr,d,lambda,tht_eb_tht_em
REAL*8 fu1,fu2

REAL*8 pr0_bar,d_bar, tht_eb_tht_el, pi


REAL*8 xg
REAL*8 u_temp(8)
REAL*8  ec





pi = 4*DATAN(1.d0)

two_sqrt2 = 2*DSQRT(2.d0)


DO i=1,n
!         U1(I)= U1(I)*DEXP( - UD*DT)! exact solution for U1
!         U2(I)= U2(I)*DEXP( - UD*DT)! exact solution for U2


  u_temp(1) = u1(i)
  u_temp(2) = u2(i)
  u_temp(3) = theta1(i)
  u_temp(4) = theta2(i)
  u_temp(5) = theta_eb(i)
  u_temp(6) = q(i)
  u_temp(7) = hs(i)
  u_temp(8) = hc(i)

  ! fu1 = - ud*u1(i)
  ! fu2 = - ud*u2(i)
  ! u1(i) = u1(i) + dt*fu1
  ! u2(i) = u2(i) + dt*fu2


!     PREDICTION




  tht_eb_tht_em = theta_eb_m_theta_em + theta_eb(i)  &
      - (two_sqrt2/pi)*(theta1(i) + alpha3*theta2(i))- q(i)
  tht_eb_tht_el = theta_eb_elbar + theta_eb(i)  &
      - (two_sqrt2/pi)*(theta1(i) + 2.d0*theta2(i))- 2.d0*q(i)

  ec=deltac*m0*tht_eb_tht_el*(two_sqrt2/pi)*(hc(i)/qr01);
  ec=deltac1*ec
  d = DMAX1(0.d0,(1.d0 + mu*(hs(i)- hc(i))/qr01) )*m0 *tht_eb_tht_em
  d=d+DMAX1(ec,0.d0)

  pr=  hd(i) + xic*hc(i) + xis*hs(i)



  ftht1 = pr - qr01 - theta1(i)/tau_r

  theta1 (i)= theta1(i) + dt* ftht1

  ftht2 = hc(i)-hs(i)-qr02 - theta2(i)/tau_r

  theta2 (i)= theta2(i) + dt* ftht2


  xg=SIN(2.d0*pi*((n/24.d0)*timel*t/(3600.d0) +i)/n- pi/2.d0)
  xg=(xg-1.d0/pi)*(pi/(pi-1.d0))
  xgtemp=(1.d0/(1.d0-pi))
  xg=DMAX1(xg,xgtemp)*s1




  fthteb =(theta_ebs_m_theta_eb+sdc*xg+thteb_st(i) -theta_eb(i))/tau_e - d/zb

  theta_eb(i)=  theta_eb(i) + dt* fthteb
  fq = d/zt - two_sqrt2*pr/pi
  q(i) = q(i) + dt*fq

  fhs = (hds(i)-hs(i))/tau_s

  hs (i)= hs(i) + dt* fhs


!   CORRECTION


  ! fu1= (fu1- ud*u1(i))/2
  ! fu2= (fu2- ud*u2(i))/2
  ! u1(i) = u_temp(1) + dt*fu1
  ! u2(i) = u_temp(2) + dt*fu2





  tht_eb_tht_em = theta_eb_m_theta_em + theta_eb(i)  &
      - (two_sqrt2/pi)*(theta1(i)+alpha3*theta2(i)) - q(i)

  tht_eb_tht_el = theta_eb_elbar + theta_eb(i)  &
      - (two_sqrt2/pi)*(theta1(i) + 2.d0*theta2(i))- 2.d0*q(i)

  ec=deltac*m0*tht_eb_tht_el*(two_sqrt2/pi)*(hc(i)/qr01);
  ec=deltac1*ec

  d = DMAX1(0.d0,(1.d0 + mu*( hs(i)-hc(i))/qr01))*m0 *tht_eb_tht_em

  d=d+DMAX1(ec,0.d0)

  pr=  hd(i) + xic*(hc(i)) + xis*( hs(i))


  ftht1 = (ftht1+pr - qr01 - theta1(i)/tau_r)/2

  theta1 (i)= u_temp(3) + dt* ftht1


  ftht2 = (ftht2 + hc(i)  - hs(i) - theta2(i)/tau_r -qr02)/2

  theta2 (i)= u_temp(4) + dt* ftht2
  xg=SIN(2.d0*pi*((n/24.d0)*(timel+dt)*t/(3600.d0) +i)/n- pi/2.d0)

  xg=(xg-1.d0/pi)*(pi/(pi-1.d0))
  xgtemp=1.d0*(1.d0/(1.d0-pi))
  xg=DMAX1(xg,xgtemp)*s1



  fthteb=(fthteb + (theta_ebs_m_theta_eb+sdc*xg+thteb_st(i)  &
      - theta_eb(i))/tau_e  - d/zb)/2

  theta_eb(i)=  u_temp(5) + dt* fthteb

  fq = (fq+d/zt - two_sqrt2*pr/pi)/2

  q(i) = u_temp(6) + dt*fq

  fhs = (fhs+(hds(i)-hs(i))/tau_s)/2
  hs(i)= u_temp(7)+dt*fhs
!         Hs(I)=Hds(i)


END DO

RETURN
END SUBROUTINE range_kuttas
end module
