module multicloud_mod
  use param_mod
  use util
  implicit none
  public  :: updatehcds, equildistr2, range_kuttas, init_multicloud, get_eqcloudfrac
  public :: T, L, alpha_bar, c, ud, theta_ebs_m_theta_eb

  private

  ! basic nondimensional scales
  real(8) :: t, c, l, alpha_bar

  integer, parameter :: nstochloc = 30

  ! parameters of mc

  real(8) nstochgl

  real(8) a,b,a0,a1,a2,a0p,a1p,a2p,alpha2,alpha3,  &
       lambdas,alpha_s,alpha_c,xis,xic,mu,ud,thetap,thetam
  real(8) tau_conv,tau_e, tau_s,tau_c,tau_r, tau_d, sdc,dtype

  real(8) tau01,tau02,tau10,tau12,tau20,tau30,tau23, r23value, times, taumult
  real(8) capebar, cape0,dry,pbar, qbar,fceq,fdeq,fseq, rstoch

  real(8) theta_eb_elbar,deltac,deltac1

  real(8) qr01,qr02,lambdabar,qcbar,dbar,  &
       theta_eb_m_theta_em,m0,theta_ebs_m_theta_eb,moist0,alpha4

  real(8) u0, cd, lcp

contains
  subroutine get_eqcloudfrac(a, b, c)
    real(8) a, b, c
    a = fceq
    b = fdeq
    c = fseq
  end subroutine get_eqcloudfrac
  subroutine init_multicloud()

    nstochgl=nstochloc
    pi=4*DATAN(1.d0)

    !     Parametrization parameters: a1=1; a2=0 ====> CAPE parametrization
    !                                 a1=0, a2=1=====> moisture parametrization


    a1=0.5D0
    a2=1.d0 -a1
    a0=2.d0

    a1p=1.d0  !Coefficient  of theta_eb in low CAPE, should be 1
    a2p=0.d0  ! Coefficient of q in low CAPE should be zero
    a0p=1.5D0
    alpha2=.1D0 ! gamma_2 in the papers (JASI,JASII,TAO,TCFD)
    alpha3=.1D0 ! alpha_2 in the papers
    alpha4=2.d0

    !     lower threshold for parameter lambda: mesuring dryness and moisteness of middle troposphere


    lambdas=0.d0

    ! Scales
    c         = dsqrt(n2) *zt / pi
    alpha_bar = zm * n2 * theta0 / g ! temperature scale
    L         = dsqrt(c / beta)
    T         = L / c


    cd=.001 ! momentum drag coefficient
    u0=2.d0 !m/s strength of turbulent fluctuations
    tau_d = 75*day! sec           Rayleigh-wind Relaxation Scale
    tau_r= 50*day !sec           Newtonian Cooling Relaxation Time Scale
    tau_s=3.d0*hour  !sec           Stratiform adjustment time Scale (not to confuse with stratiform time scale which is given by tau_conv
    tau_conv=  2*hour !sec           Convective time scale
    tau_c=3*hour   !sec   Congestus adjustment time scale (not to confuse with congestus time scale which is given by tau_conv/alpha_c )



    ! RCE fixed by QR01,theta_ebs_m_theta_eb,theta_eb_m_theta_em,thetap,thetam




    qr01 = 1.d0/day!/2/dsqrt(2.d0)       ! prescribed radiative cooling associated with first baroclinic

    alpha_s=0.25D0            ! ratio of deep convective and stratiform heating rates

    alpha_c=0.1D0  !ratio of deep convective and congestus time scales

    mu= 0.25D0   ! contribution of lower  tropospheric cooling (H_s - H_c) to downdrafts
    xic=0.0D0
    xis=0.4D0


    deltac=(16-3*pi*xic)/(16+3*pi*xic)

    theta_ebs_m_theta_eb = 10.d0 !K discrepancy between boundary layer theta_e and its saturation value

    tau_e= (theta_ebs_m_theta_eb/qr01)*(zb/zt)*pi/2.d0/DSQRT(2.d0)! Evaporative time scale

    theta_eb_m_theta_em=11.d0!K discrepancy between boundary layer and middle tropospheric theta_e''s
    theta_eb_elbar=0.d0!K discrepancy between boundary layer and lower middle tropospheric theta_e''s

    deltac1=0.d0




    thetap=20.d0!
    thetam=10.d0 !K thresholds moistening and drying of middle troposphere

    a=(1.d0-lambdas)/(thetap-thetam)
    b=lambdas-a*thetam            !linear fit coefficients of LAMBDA


    ! Stochastic Params

    CAPE0 = 400 ! J / Kg
    MOIST0 = 30 ! K
    rstoch = 2 * zp* cp * gammam / theta0 *alpha_bar / c /c



    IF(theta_eb_m_theta_em <= thetam) THEN
       lambdabar=lambdas
    ELSE IF(theta_eb_m_theta_em < thetap) THEN
       lambdabar=a*theta_eb_m_theta_em+b
    ELSE
       lambdabar=1.d0
    END IF



    !    stochastic RCE parameters
    ! TODO : refactor m0, rstoch, alpha_Bar, l, c into calculate_rce
    
    taumult = 1.0
    ! FMK13 Taus
    tau01 = 1.0d0
    tau02 = 3.0d0
    tau10 = 1.0d0
    tau12 = 1.0d0
    tau20 = 3.0d0
    tau30 = 5.0d0
    tau23 = 3.0d0


    tau01=tau01*taumult
    tau02=tau02*taumult
    tau10=tau10*taumult
    tau12=tau12*taumult
    tau20=tau20*taumult
    tau30=tau30*taumult
    tau23= tau23*taumult
    call nondimensionalize_params
    call calculate_rce()

  end subroutine init_multicloud

  subroutine nondimensionalize_params()



    t=l/c            !~ 8.33 hours, time scale


    lcp= (1000/700)**(287.4/1004)*(2.504E6)/1004/alpha_bar

    qcbar=qbar

    a =a *alpha_bar

    theta_eb_m_theta_em=theta_eb_m_theta_em/alpha_bar
    theta_ebs_m_theta_eb =  theta_ebs_m_theta_eb/alpha_bar
    theta_eb_elbar =  theta_eb_elbar/alpha_bar


    thetap=thetap/alpha_bar
    thetam=thetam/alpha_bar

    qr01 = qr01*t/alpha_bar

    !       QR02 = QR02*T/ALPHA_BAR


    ud = cd*u0*t/zb + t/tau_d

    zb=zb/l
    zt=zt/l
    zm=zm/l

    tau_r=tau_r/t
    tau_s=tau_s/t
    tau_c=tau_c/t
    tau_conv = tau_conv/t


    ! Stochastic Params

    CAPE0 = CAPE0 / c / c
    MOIST0 = MOIST0 / alpha_bar

  end subroutine nondimensionalize_params

  subroutine calculate_rce()
    ! This code doesn't handle congestus detrainment
    real *8 :: hsbar, hcbar, hdbar, tspi
    call calculate_srce()

    hsbar = fseq * alpha_s * dsqrt(capebar) /zm
    hcbar = fceq * alpha_c * dsqrt(capebar) /zm
    hdbar = fdeq * dsqrt(capebar) / zm


    qr02 = Hcbar - Hsbar

    pbar = Hdbar + Hsbar * xis + Hcbar * xic

    ! This dbar is different from ColumnStandard.m, but correct I think
    dbar = 2.d0 * dsqrt(2.0d0)/pi * pbar * zt

    m0   = (1.0d0 + mu * ( Hsbar - Hcbar) / qr01 )* theta_eb_m_theta_em
    m0   = dbar /m0

    tau_e = theta_ebs_m_theta_eb * zb / dbar

  end subroutine calculate_rce


  subroutine check_srce()
    real *8 x(4), fvec(4)

    x = (/capebar, fceq, fdeq, fseq/)
    call fcn(4, x, fvec, 0)

    ! print *, 'Should be close to zero'
    ! print *,fvec
    ! print *
  end subroutine check_srce

  subroutine calculate_srce()

    integer info, i

    real *8 fvec(4), x(4)
    real *8 wa(200)
    real *8 :: tol = 1d-16

    x = (/(.1d0, i = 1,4)/)
    call hybrd1(fcn,4,x,fvec,tol,info,wa,200)

    capebar = x(1)
    fceq    = x(2)
    fdeq    = x(3)
    fseq    = x(4)
  end subroutine calculate_srce

  subroutine fcn(n,x,fvec,iflag)
    integer n,iflag
    double precision x(n),fvec(n)
    real *8 cd,cl,d,r01,r02,r10,r20,r23,r12,r30
    real *8 sig_clear

    capebar = x(1)
    qbar =  dsqrt( capebar )/ zm

    fceq    = x(2)
    fdeq    = x(3)
    fseq    = x(4)

    sig_clear = 1.0d0 - fceq - fdeq - fseq

    cd = capebar /cape0
    cl = capebar / cape0
    d  = theta_eb_m_theta_em / moist0

    call equildistr2(cd,cl,d,r01,r02,r10,r20,r23,r12,r30)

    fvec(1) = sig_clear * r01 - fceq * (r10 + r12)
    fvec(2) = sig_clear * r02 + fceq * r12 - fdeq * (r20 + r23)
    fvec(3) = fdeq * r23 - fseq * r30
    fvec(4) = QR01 - qbar * (fdeq + fseq * xis * alpha_s  &
         + fceq * xic * alpha_c)

  end subroutine fcn


  SUBROUTINE range_kuttas(u1,u2,theta1,theta2,theta_eb,q,hs,hc,hd,  &
       n,dt, thteb_st,timel,hds, moiststab)
    use param_mod
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
    REAL*8, INTENT(OUT)                      :: moiststab(n)
    INTEGER :: i

    REAL*8 ftht1,ftht2,fthteb,fq,fhs,fhc,two_sqrt2,xgtemp
    REAL*8 pr0,pr1,pr,d,lambda,tht_eb_tht_em
    REAL*8 fu1,fu2

    REAL*8 pr0_bar,d_bar, tht_eb_tht_el


    REAL*8 xg
    REAL*8 u_temp(8)
    REAL*8  ec



    two_sqrt2 = 2*DSQRT(2.d0)


    DO i=1,n
       ! U1(I)= U1(I)*DEXP( - UD*DT)! exact solution for U1
       ! U2(I)= U2(I)*DEXP( - UD*DT)! exact solution for U2


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

       moiststab(i) = -tht_eb_tht_em

       ec=deltac*m0*tht_eb_tht_el*(two_sqrt2/pi)*(hc(i)/qr01);
       ec=deltac1*ec
       d = DMAX1(0.d0,(1.d0 + mu*(hs(i)- hc(i))/qr01) )*m0 *tht_eb_tht_em
       d=d+DMAX1(ec,0.d0)

       pr=  hd(i) + xic*hc(i) + xis*hs(i)



       ftht1 = pr - qr01 - theta1(i)/tau_r

       theta1 (i)= theta1(i) + dt* ftht1

       ftht2 = hc(i)-hs(i)-qr02 - theta2(i)/tau_r

       theta2 (i)= theta2(i) + dt* ftht2





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

       fthteb=(fthteb + (theta_ebs_m_theta_eb+thteb_st(i)  &
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

  SUBROUTINE equildistr2(cd,cl,d,r01,r02,r10,r20,r23,r12,r30)
    ! Calculate the transition rates given the dryness, lower level
    ! CAPE and middle level CAPE
    use param_mod
    IMPLICIT NONE
    REAL*8, INTENT(IN OUT)                   :: cd
    REAL*8, INTENT(IN OUT)                   :: cl
    REAL*8, INTENT(IN OUT)                   :: d
    REAL*8, INTENT(OUT)                      :: r01
    REAL*8, INTENT(OUT)                      :: r02
    REAL*8, INTENT(OUT)                      :: r10
    REAL*8, INTENT(OUT)                      :: r20
    REAL*8, INTENT(OUT)                      :: r23
    REAL*8, INTENT(OUT)                      :: r12
    REAL*8, INTENT(OUT)                      :: r30
    REAL*8  dryness,il,tht_eb_tht_emloc




    r01 = gammabb(cl)*gammabb(d)/tau01
    r02 = gammabb(cd)*(1.d0-gammabb(d))/tau02
    r10 = gammabb(d)/tau10
    r12 = gammabb(cd)*(1.d0-gammabb(d))/tau12
    r20 = (1.d0-gammabb(cd))/tau20
    r30 = 1.0/tau30
    r23 = 1.d0/tau23
  END SUBROUTINE equildistr2



  !   service subroutine for Birthd. funct

  SUBROUTINE findindex(rmat,rsum,test,rindex)
    ! This subroutine chooses the appropriate transition
    ! given a uniform random number test
    !
    ! Arguments:
    !
    ! rmat is the concatenated rates
    ! rsum is a dummy variable
    ! rindexis the chosen transition

    IMPLICIT NONE
    REAL*8, INTENT(IN)                       :: rmat(7)
    REAL*8, INTENT(OUT)                      :: rsum(7)
    REAL*8, INTENT(IN)                       :: test
    INTEGER, INTENT(OUT)                     :: rindex

    REAL*8  rsumtemp(7)
    INTEGER :: temp1,il

    ! Take the cumsum of rmat
    rsum(1) = rmat(1)
    DO  il = 2, 7
       rsum(il) = rsum(il-1) + rmat(il)
    END DO

    ! Shift cumsumed rsum by test
    rsumtemp=rsum/rsum(7)-[1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0]*test


    ! rindex is the interval surrounding zero
    il=1
    temp1=0.d0


    DO  il = 1, 7
       IF(temp1 == 0.d0) THEN
          IF( rsumtemp(il) > 0.d0) THEN
             rindex=il
             temp1=1
          END IF
       END IF
    END DO

  END SUBROUTINE findindex

  !     birth death routine
  SUBROUTINE   birthdeath2(fcloc,fdloc,fsloc,cd,cl,d,dt, nstoch)

    REAL*8, INTENT(OUT)                      :: fcloc
    REAL*8, INTENT(OUT)                      :: fdloc
    REAL*8, INTENT(OUT)                      :: fsloc
    REAL*8, INTENT(IN OUT)                   :: cd
    REAL*8, INTENT(IN OUT)                   :: cl
    REAL*8, INTENT(IN OUT)                   :: d
    REAL*8, INTENT(IN)                       :: dt
    REAL*8, INTENT(IN)                       :: nstoch

    REAL*8 useed,ubits,dunib,dustar,duni
    REAL*8 r01,r02,r10,r20,r23,r12,r30
    REAL*8 time,count,t,tloc
    REAL*8  lr01,lr02,lr10,lr20,lr23,lr12,lr30
    REAL*8 rsum(7), rmat(7),test
    REAL*8  diff(3), n
    REAL*8 tlocal,timeloc,nsite, temp3
    INTEGER :: rindex, temp2,ind,clearsky,cloud(3),diffv(7,3)

    rindex=0

    nsite = nstoch*nstoch

    !     Conditional rates
    CALL equildistr2(cd,cl,d,r01,r02,r10,r20,r23,r12,r30)
    !     Define cloud state

    cloud(1) = DNINT(nsite*fcloc)
    cloud(2) = DNINT(nsite*fdloc)
    cloud(3) = DNINT(nsite*fsloc)


    !     define the seven possible transitions
    diffv(1,1:3) = [1,0,0]
    diffv(2,1:3) = [-1,0,0]
    diffv(3,1:3) = [-1,1,0]
    diffv(4,1:3) = [0,1,0]
    diffv(5,1:3) = [0,-1,0]
    diffv(6,1:3) = [0,-1,1]
    diffv(7,1:3) = [0,0,-1]

    !     We will use Gilespie's algorithm and iterate until we have exceeded the
    !     time step, DT
    timeloc = 0
    count = 0
    clearsky = DNINT(nsite) - sum(cloud)

50  IF (timeloc < dt) THEN
       !      write(*,*) timeloc, dt
       count = count+1
       !   Calculate the number of clearsky elements
       clearsky = nsite - sum(cloud)
       !     Absolute rates

       lr01 = r01*clearsky
       lr10 = r10*cloud(1)
       lr12 = r12*cloud(1)
       lr02 = r02*clearsky
       lr20 = r20*cloud(2)
       lr23 = r23*cloud(2)
       lr30 = r30*cloud(3)

       rmat = [lr01,  lr10,  lr12,  lr02,  lr20,  lr23,  lr30]

       test=duni()

       CALL  findindex(rmat,rsum,test,rindex)

       IF(rsum(7) < 0.d0) THEN
          WRITE(*,*) "RSUM ERROR"
          stop
       END IF

       IF(rsum(7) == 0.d0) THEN
          tloc = 2.d0*dt
       END IF

       IF(rsum(7) > 0.d0) THEN
          tloc = -1.d0*LOG(duni())/rsum(7)
       END IF





       !     Calculate the time until the next transition



       timeloc = timeloc+tloc

       !      reject transition if time exceeds DT
       IF(timeloc <= dt) THEN

          !      Which transition occurs?

          ind = rindex



          !    write(*,*) rindex

          diff = diffv(ind,:)




          !       New cloud array

          cloud = cloud + diff
       END IF


       GO TO 50
    END IF
    fcloc = cloud(1)/nsite
    fdloc = cloud(2)/nsite
    fsloc = cloud(3)/nsite


    RETURN
  END SUBROUTINE   birthdeath2

  SUBROUTINE updatehcds(fcls,fdls,fsls, u1, u2, theta1,theta2,theta_eb,q, hds,hc,hd,n, dx, dt)
    use param_mod

    REAL*8, INTENT(IN OUT)                   :: fcls(n)
    REAL*8, INTENT(IN OUT)                   :: fdls(n)
    REAL*8, INTENT(IN OUT)                   :: fsls(n)
    REAL*8, INTENT(IN OUT)                   :: theta1(n)
    REAL*8, INTENT(IN)                       :: theta2(n)
    REAL*8, INTENT(IN)                       :: theta_eb(n)
    REAL*8, INTENT(IN)                       :: q(n)
    REAL*8, INTENT(OUT)                      :: hds(n)
    REAL*8, INTENT(OUT)                      :: hc(n)
    REAL*8, INTENT(OUT)                      :: hd(n)
    INTEGER, INTENT(IN)                      :: n
    REAL*8, INTENT(IN)                   :: dt, dx
    INTEGER :: i,iL
    REAL*8 u1(n),u2(n), w1(n), w2(n)
    REAL*8 ftht1,ftht2,fthteb,fq,fhs,fhc,two_sqrt2,xgtemp

    real(8) :: hmloc, cd, cl, tht_eb_tht_emloc, fcloc, fdloc, fsloc, dryness



    two_sqrt2 = 2*DSQRT(2.d0)
    hmloc=5.d0*zt/16.d0
    !      Rstoch=0.0d0


    DO  il = 1, n

       cd= capebar
       cd=cd+rstoch*(theta_eb(il)-1.7D0*(theta1(il)+alpha2*theta2(il)))
       cd = cd/cape0
       cl= capebar
       cl=cl+rstoch*(theta_eb(il)-1.7D0*(theta1(il)+alpha4*theta2(il)))
       cl = cl / cape0

       cd=DMAX1(cd,0.d0)
       cl=DMAX1(cl,0.d0)

       tht_eb_tht_emloc = theta_eb_m_theta_em + theta_eb(il)  &
            - (two_sqrt2/pi)*(theta1(il) + alpha3*theta2(il))- q(il)

       fcloc=fcls(il)
       fdloc=fdls(il)
       fsloc=fsls(il)


       dryness=tht_eb_tht_emloc/moist0
       CALL birthdeath2(fcloc,fdloc,fsloc, cd ,cl, dryness,dt, nstochgl)

       if (isnan(fcloc)) then
          print *, 'Error: Birthdeath yields NaN...stopping'
          stop
       endif

       hd(il)=a1*theta_eb(il)+a2*q(il)
       hd(il)=hd(il)-a0*(theta1(il)+alpha2*theta2(il))
       hd(il)=(1.d0/(tau_conv))*hd(il)
       hd(il)=hd(il)+fdeq*DSQRT(capebar)/hmloc
       hd(il)=DMAX1(hd(il),0.d0)*(fdloc/fdeq)


       hds(il)=alpha_s*hd(il)*(fsloc/fseq)
       hc(il) =  fcloc*alpha_c*SQRT(DMAX1(cl,0.d0))/hmloc


       fcls(il) =fcloc
       fdls(il) =fdloc
       fsls(il) =fsloc




    END DO



  END SUBROUTINE updatehcds

  function gammabb(xtemp)
    real *8 xtemp
    real *8 gammabb
    gammabb = 1.d0 - exp (-dmax1(xtemp, 0.0d0))
  end function gammabb

end module multicloud_mod
