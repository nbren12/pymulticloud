module stochastic 
  use util
  private

  public  :: updatehcds, w0, equildistr2

contains



  !   equilibrium distribution

  SUBROUTINE equildistr2(cd,cl,d,r01,r02,r10,r20,r23,r12,r30)
    ! Calculate the transition rates given the dryness, lower level
    ! CAPE and middle level CAPE
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
    REAL*8 tau01,tau02,tau10,tau12,tau20,tau30,tau23, r23value, times
    REAL*8 capebar,cape0,dry,pbar,qbar,fceq,fdeq,fseq,rstoch,nstochgl
    REAL*8  dryness,il,tht_eb_tht_emloc,moist0


    COMMON /tauscales/ capebar, cape0,dry,pbar, qbar,fceq,fdeq,fseq
    COMMON /tauscales/  tau01,tau02,tau10,tau12, rstoch
    COMMON /tauscales/  tau20,tau30,tau23, r23value, times,nstochgl



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


    REAL*8 pr0,pr1,pr,d,lambda,tht_eb_tht_em,alpha3,hsbar, hcbar
    REAL*8 a0p,a1p,a2p,cconv,alpha4,theta_ebs_m_theta_eb

    REAL*8 a,b,a0,a1,a2,alpha2,lambdas,alpha_s,alpha_c,mu,ud
    REAL*8 thetap,thetam,zb,zt,fu1,fu2,xis,xic

    REAL*8 qr01,qr02,lambda_bar,pr0_bar,d_bar, theta_eb_m_theta_em,m0,moist0



    REAL*8 tau_conv,tau_e, tau_s,tau_c,tau_r,pi,timel,t,xg, sdc,dtype
    REAL*8 u_temp(8)
    REAL*8 thteb_st(n),nstochgl

    REAL*8 tau01,tau02,tau10,tau12,tau20,tau30,tau23, r23value, times
    REAL*8 capebar, cape0,dry,pbar, qbar,fceq,fdeq,fseq, rstoch
    REAL*8 cd, cl, wd, wl, dryness,tht_eb_tht_emloc


    REAL*8 fcloc, fdloc, fsloc, hmloc,lambdabar,qcbar,dbar

    real*8 conv1, conv2
    COMMON/const_par/a,b,a0,a1,a2,a0p,a1p,a2p,alpha2,alpha3,  &
         lambdas,alpha_s,alpha_c,xis,xic,mu,ud,thetap,thetam,zb,zt

    COMMON /char_times/tau_conv,tau_e, tau_s,tau_c,tau_r,t, sdc,dtype
    COMMON /tauscales/ capebar, cape0,dry,pbar, qbar,fceq,fdeq,fseq
    COMMON /tauscales/  tau01,tau02,tau10,tau12, rstoch
    COMMON /tauscales/  tau20,tau30,tau23, r23value, times,nstochgl
    COMMON /rce_values/qr01,qr02,lambdabar,qcbar,dbar,  &
         theta_eb_m_theta_em,m0,theta_ebs_m_theta_eb,moist0,alpha4


    two_sqrt2 = 2*DSQRT(2.d0)
    pi=4*DATAN(1.d0)
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

  subroutine fmk13_C(cl, cd, theta1, theta2, theta_eb)
    real*8 cl, cd, theta1, theta2, theta_eb

    cd= capebar
    cd=cd+rstoch*(theta_eb-1.7D0*(theta1+alpha2*theta2))
    cd = cd/cape0
    cl= capebar
    cl=cl+rstoch*(theta_eb-1.7D0*(theta1+alpha4*theta2))
    cl = cl / cape0
  end subroutine fmk13_C

end module stochastic


