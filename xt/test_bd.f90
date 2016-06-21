
! Code converted using TO_F90 by Alan Miller
! Date: 2014-04-30  Time: 02:18:49

PROGRAM test_birth_death
IMPLICIT NONE
INTEGER :: n,i, niter, iter,k,nout, iout,inrg
REAL*8 taumult ,nstochgl,nstochloc
PARAMETER(n=1040)
PARAMETER(nstochloc=30.d0)
PARAMETER(taumult=1.d0)

!     PROGNOSTIC VARIABLES
REAL*8 u1(n),u2(n),theta1(n),theta2(n),  &
    q(n),theta_eb(n),hd(n),hs(n),hc(n),hds(n)
REAL*8 uc(5,-1:n+2)
REAL*8 thteb_st(n), asst,lsst, sdc,dtype






REAL*8 q_tld,alpha_tld,lmd_tld,tempg,dctime,dcint,s1

PARAMETER(alpha_tld=0.1D0,lmd_tld=0.8D0,q_tld=0.9D0)

REAL*8 dt,dx,p, umax,time,dt_max,twave_count2
REAL*8 alpha_bar,c,l,t,theta_d,tend,tout,tenergy





! Parameters for precipitation relaxation scheme

REAL*8 a1,a2,a0,alpha2,fs,fd,fc,deltas,deltac,alpha4
!      PARAMETER(a1=.9,a2=.1,a0=1.7,alpha2=1.d0)


REAL*8 hour, minute, day, km, tempouttime
PARAMETER(hour=3600.d0, minute=60.0D0, day=86400.0D0, km=1000.d0)

REAL*8 pi,n2,gamma,cp,omega,r, theta0,g, beta, lcp,EQ
REAL*8 zt,zb,u0, cd, tau_d, tau_r, tau_s, tau_conv, tau_c, lambdas
REAL*8 qr01, qr02, alpha_s, alpha_c, xic, xis, mu
REAL*8 theta_ebs_m_theta_eb, tau_e,theta_eb_m_theta_em,theta_embar

REAL*8 thetap, thetam, a, b, m0, ud ,theta_ebs, theta_ebbar

REAL*8   qcbar, dbar,  lambdabar,alpha3,a0p,a1p,a2p
REAL*8 u1_av(n),u2_av(n),tht1_av(n),tht2_av(n),tht_eb_av(n),  &
    q_av(n),hs_av(n),hc_av(n),twave_count,twave,twave_out, pr0_av(n),hd_av(n),  &
    lambda_av(n),lambda,tht_eb_tht_em,pr0,moist0

!     stochastic variable code
REAL*8 tau01,tau02,tau10,tau12,tau20,tau30,tau23, r23value, times
REAL*8 capebar, cape0,dry,pbar, qbar,fceq,fdeq,fseq, rstoch

REAL*8 fcls(n),fdls(n),fsls(n),theta_eb_elbar,deltac1


COMMON/const_par/a,b,a0,a1,a2,a0p,a1p,a2p,alpha2,alpha3,  &
    lambdas,alpha_s,alpha_c,xis,xic,mu,ud,thetap,thetam,zb,zt

COMMON /rce_values/qr01,qr02,lambdabar,qcbar,dbar,  &
    theta_eb_m_theta_em,m0,theta_ebs_m_theta_eb,moist0,alpha4

COMMON /char_times/tau_conv,tau_e, tau_s,tau_c,tau_r,t, sdc,dtype

COMMON /tauscales/ capebar, cape0,dry,pbar, qbar,fceq,fdeq,fseq
COMMON /tauscales/  tau01,tau02,tau10,tau12, rstoch
COMMON /tauscales/  tau20,tau30,tau23, r23value, times,nstochgl


COMMON /congd/ theta_eb_elbar,deltac,deltac1

NAMELIST /deterministic/ asst, lsst, theta_eb_m_theta_em
NAMELIST /tauscales/ capebar, cape0,dry,pbar, qbar,fceq,fdeq,fseq
NAMELIST /tauscales/  tau01,tau02,tau10,tau12, rstoch
NAMELIST /tauscales/  tau20,tau30,tau23, r23value, times,nstochgl

!     on fly averaging

REAL*8 u1a(n),u2a(n),theta1a(n),theta2a(n),theta_eba(n)
REAL*8 qa(n),hsa(n),hca(n),hda(n),fclsa(n),fdlsa(n),fslsa(n)
REAL*8 twave_outa,twave_count2a,tenergya,twa


!     on fly averaging



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

!     Benchmark of standard parameters



cp=1000.d0     !J/kg/K
gamma=1.7D0   !none

omega= 2*pi/(24*hour)! the angular velocity of the earth and
r=6378*km           ! radius of the earth
theta0=300.d0 !K background  reference temperature
g=9.8D0 !m/s gravitational acceleration
beta=2*omega/r      !1/s/m
EQ=40000*km!  Earth's peremeter at the Equator
n2=.0001D0         !Vaissala buoyancy frequency squared


zt=16.00*km !topospheric height
zb=500.d0 !meters boundary layer depth
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
fs = 0.4D0 ! stratiform fraction  of total  surface rain fall.
fd = 0.6D0 !   deep convective fraction of total  surface rain fall.
fc = 0.d0 ! congestus fraction of total   surface precipitation


xic= fc/fd/alpha_c
xis = fs/fd/alpha_s
xic=0.0D0
xis=0.4D0


PRINT*,'xic=',xic
PRINT*,'xis=',xis


deltas=(16-3*pi*xis)/(16+3*pi*xis)

deltac=(16-3*pi*xic)/(16+3*pi*xic)

PRINT*,'deltas=',deltas,'deltac=',deltac

111  CONTINUE

!       print*, 'Enter 0 to shut mu off 1 otherwise:'

!       read*, mu



theta_ebs_m_theta_eb = 10.d0 !K discrepancy between boundary layer theta_e and its saturation value

tau_e= (theta_ebs_m_theta_eb/qr01)*(zb/zt)*pi/2.d0/DSQRT(2.d0)! Evaporative time scale

theta_eb_m_theta_em=11.d0!K discrepancy between boundary layer and middle tropospheric theta_e''s
theta_eb_elbar=0.d0!K discrepancy between boundary layer and lower middle tropospheric theta_e''s

deltac1=0.d0




thetap=20.d0!
thetam=10.d0 !K thresholds moistening and drying of middle troposphere

a=(1.d0-lambdas)/(thetap-thetam)
b=lambdas-a*thetam            !linear fit coefficients of LAMBDA


OPEN(69,FILE="SRCE.txt",STATUS='old')

PRINT*,'Evaporative time tau_e=', tau_e/3600/24,' days'
IF(theta_eb_m_theta_em <= thetam) THEN
  lambdabar=lambdas
  
ELSE IF(theta_eb_m_theta_em < thetap) THEN
  lambdabar=a*theta_eb_m_theta_em+b
  
ELSE
  lambdabar=1.d0
END IF







!  Reference scales

 ! m isture: ratio of latent heat of vaporization and heat capacity
!times (p/p_0)^{-.28), p=700 hPa, P_0=1000 normalized by alpha_bar unit of temperature.


!    NON DIMENTIONALIZATION OF RCE VALUES AND PARAMETERS USED BEYOND THIS


!    stochastic RCE parameters
WRITE (*,*) 'rcefile'
READ(69,*) qbar
READ(69,*) dbar
READ(69,*) pbar
READ(69,*) qr02
READ(69,*) tau_e
READ(69,*) m0
READ(69,*) cape0
READ(69,*) capebar
READ(69,*) fceq
READ(69,*) fdeq
READ(69,*) fseq
READ(69,*) rstoch


READ(69,*) tau01
READ(69,*) tau02
READ(69,*) tau10
READ(69,*) tau12
READ(69,*) tau20
READ(69,*) tau30
READ(69,*) tau23

READ(69,*) alpha_bar
READ(69,*) l
READ(69,*) c
READ(69,*) moist0


open(8, file="input.nml", status='OLD')
read(8, nml = deterministic)
write (*, NML=deterministic)
close(8)

tau01=tau01*taumult
tau02=tau02*taumult
tau10=tau10*taumult
tau12=tau12*taumult
tau20=tau20*taumult
tau30=tau30*taumult
tau23= tau23*taumult




t=l/c            !~ 8.33 hours, time scale
theta_d=alpha_bar !~15 degree K; temperature scale


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
WRITE (*,*) 'zb'
WRITE (*,*) zb

tau_r=tau_r/t
tau_s=tau_s/t
tau_c=tau_c/t
tau_conv = tau_conv/t










DO i=1,n
  fcls(i)=fceq;
  fdls(i)=fdeq;
  fsls(i)=fseq;
END DO


WRITE(*,*) 'tau10'
WRITE(*,*) tau10
PRINT*,'m0',m0, 'mu/Qc',mu
!        stop

!     COMMON/CONST_PAR/A,B,A0,A1,A2,ALPHA2,LAMBDAS,ALPHA_S,ALPHA_C,MU,
!     x     UD,thetap,thetam,ZB,ZT
!      COMMON /RCE_VALUES/QR01,QR02,LAMBDABAR,QcBAR,DBAR
!     X     THETA_EB_M_THETA_EM,M0
!      COMMON /CHAR_TIMES/TAU_CONV,TAU_E, TAU_S,TAU_C,tau_r



!    GRID and MAXIMUM TIME STEP


!       print*,'ENTER WAVENUMBER k'
!       read*, k
!        P=EQ/L


dt_max =2.d0*minute/t
dt_max =1.d0*minute/t



tend =  500.00D0*day/t !day/T

niter=tend/dt_max
nout=1

tout=tend/nout

tenergya= dt_max*3.d0
twave_outa=tend-100.d0*day/t

tenergy= (3.d0)*hour/t
twave_out=400.d0*day/t



CALL initial_data(u1,u2,theta1,theta2,theta_eb,q,hs,hc,n,dx,l,p)
!(u1,u2,theta1,theta2,theta_eb,q,hs,hc,N,DX,L,P)



!        Prescribed STT with non-zero gradient

!        ASST : strenght of warm pool heating

!       LSST : width of warm pool.

asst=0.d0;
lsst=EQ/l/8;
!      Strength of the diurnal cycle in K
dctime=0.d0*day/t;
dcint=0.01D0*day/t;
sdc=0.0D0/theta_d;
dtype=2;








! CFL
umax=1.d0
DO i=1,n
  umax = DMAX1(umax, (DABS(q(i)+q_tld)+DABS( alpha_tld*q(i) +  &
      lmd_tld*q_tld)+ DABS(u1(i)+alpha_tld*u2(i)))/5.d0)
END DO




dt = DMIN1(0.9*dx/umax,dt_max)


DO i=1,n
  uc(1,i) = u1(i)
  uc(2,i) = u2(i)
  uc(3,i) = theta1(i)
  uc(4,i) = theta2(i)
  uc(5,i) = q(i)
END DO



tempg=dt
do i=1,1000
CALL updatehcds(fcls,fdls,fsls,theta1,theta2,theta_eb,q,hds,hc,hd  &
    ,n,2.d0*dt*t/(hour) )
enddo
dt=tempg






END PROGRAM test_birth_death


!     MAIN Program Ends here


SUBROUTINE initial_data(u1,u2,theta1,theta2,theta_eb,q,hs,hc, n,dx,l,p)
IMPLICIT NONE

REAL*8, INTENT(OUT)                      :: u1(n)
REAL*8, INTENT(OUT)                      :: u2(n)
REAL*8, INTENT(OUT)                      :: theta1(n)
REAL*8, INTENT(OUT)                      :: theta2(n)
REAL*8, INTENT(OUT)                      :: theta_eb(n)
REAL*8, INTENT(OUT)                      :: q(n)
REAL*8, INTENT(OUT)                      :: hs(n)
REAL*8, INTENT(OUT)                      :: hc(n)
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(OUT)                      :: dx
REAL*8, INTENT(IN)                       :: l
REAL*8, INTENT(OUT)                      :: p
INTEGER :: i, k, kn,idata,j,n1
REAL*8  pi, x
CHARACTER (LEN=10) :: file_name

REAL*8 u1re,u1im, u2re,u2im, tht1re,tht1im,tht2re,tht2im
REAL*8 tht_ebre,tht_ebim,qre,qim,hsre,hsim,hcre,hcim
REAL*8 epsilonz
REAL*8, PARAMETER :: epsilon=.01D0

pi=4*DATAN(1.d0)


!     TYPE OF INITIAL DATA


PRINT*,'ENTER TYPE OF INITIAL DATA: 1: VALIDATION,  &
    2:RANDOM PERTURBATION, 3: to read form a file'
!READ*, idata
idata = 1
SELECT CASE ( idata )
  CASE (    1)
    GO TO 10
  CASE (    2)
    GO TO 20
  CASE (    3)
    GO TO 30
END SELECT
10   CONTINUE
PRINT*,'ENTER NAME OF FILE CONTAINING MODE DATA'
!     READ*, FILE_NAME
file_name='deep_moist'
OPEN(UNIT=21,FILE=file_name,STATUS='OLD')


READ(21,*) k !WAVENUMBER
READ(21,*) u1re,u1im
READ(21,*) u2re,u2im
READ(21,*) tht1re,tht1im
READ(21,*) tht2re,tht2im
READ(21,*) tht_ebre,tht_ebim
READ(21,*) qre,qim
READ(21,*) hsre,hsim
READ(21,*) hcre,hcim
!      DATA U1RE,U1IM/1.387778780781446d-16,-4.367406655446615d-01/
!      DATA U2RE,U2IM/-3.122502256758253e-17 , -8.552912998553463e-03/
!      DATA THT1RE,THT1IM/8.171778296182428e-01, 0.000000000000000e+00/
!      DATA  THT2RE,THT2IM/1.600320609567062e-02, 8.673617379884035e-18/
!      DATA  THT_EBRE,THT_EBIM/1.063447982491118e-01,
!     x     -5.551115123125783e-17/
!      DATA  QRE,QIM/ -1.984886156582951e-01 ,-1.387778780781446e-17/
!      DATA  HSRE,HSIM/3.007354057195986e-01 ,8.326672684688674e-17/
!      DATA  HCRE,HCIM/0.000000000000000e+00, 0.000000000000000e+00/


CLOSE(21)
p= 40000000.d0/l/k
dx=p/n
DO i=1,n
  x= (i-1)*dx*2*pi/p
  u1(i) = epsilon*(u1re*DCOS(x) - u1im*DSIN(x))
  u2(i) = epsilon*(u2re*DCOS(x) - u2im*DSIN(x))
  theta1(i) = epsilon*(tht1re*DCOS(x) - tht1im*DSIN(x))
  theta2(i) = epsilon*(tht2re*DCOS(x) - tht2im*DSIN(x))
  theta_eb(i) = epsilon*(tht_ebre*DCOS(x) - tht_ebim*DSIN(x))
  q(i) = epsilon*(qre*DCOS(x) - qim*DSIN(x))
  hs(i) = epsilon*(hsre*DCOS(x) - hsim*DSIN(x))
  hc(i) = epsilon*(hcre*DCOS(x) - hcim*DSIN(x))
  
END DO
GO TO 40
20   CONTINUE





PRINT*,'ENTER NAME OF FILE CONTAINING MODE DATA'
!     READ*, FILE_NAME
file_name='deep_moist'
OPEN(UNIT=21,FILE=file_name,STATUS='OLD')

READ(21,*) k !WAVENUMBER
READ(21,*) u1re,u1im
READ(21,*) u2re,u2im
READ(21,*) tht1re,tht1im
READ(21,*) tht2re,tht2im
READ(21,*) tht_ebre,tht_ebim
READ(21,*) qre,qim
READ(21,*) hsre,hsim
READ(21,*) hcre,hcim

!      DATA U1RE,U1IM/1.387778780781446d-16,-4.367406655446615d-01/
!      DATA U2RE,U2IM/-3.122502256758253e-17 , -8.552912998553463e-03/
!      DATA THT1RE,THT1IM/8.171778296182428e-01, 0.000000000000000e+00/
!      DATA  THT2RE,THT2IM/1.600320609567062e-02, 8.673617379884035e-18/
!      DATA  THT_EBRE,THT_EBIM/1.063447982491118e-01,
!     x     -5.551115123125783e-17/
!      DATA  QRE,QIM/ -1.984886156582951e-01 ,-1.387778780781446e-17/
!      DATA  HSRE,HSIM/3.007354057195986e-01 ,8.326672684688674e-17/
!      DATA  HCRE,HCIM/0.000000000000000e+00, 0.000000000000000e+00/


CLOSE(21)
PRINT*,'Enter 1 for modulated sine wave packet'
PRINT*,'Enter 2 for localized--single wave'

!      read*, idata
idata=2



SELECT CASE ( idata )
  CASE (    1)
    GO TO 21
  CASE (    2)
    GO TO 22
END SELECT
21   CONTINUE
kn=k
n1=n/kn


p= 40000000.d0/l/kn
dx=p/n1


DO i=1,kn
  x= (i-1)*dx*2*pi/p
  epsilonz = epsilon*DSIN(x)
  
  DO j=1,n1
    
    x= ((i-1)*n1+j-1)*dx*2*pi/p
    k=n1*(i-1)+j
    u1(k) = epsilonz*(u1re*DCOS(x) - u1im*DSIN(x))
    u2(k) = epsilonz*(u2re*DCOS(x) - u2im*DSIN(x))
    theta1(k) = epsilonz*(tht1re*DCOS(x) - tht1im*DSIN(x))
    theta2(k) = epsilonz*(tht2re*DCOS(x) - tht2im*DSIN(x))
    theta_eb(k) = epsilonz*(tht_ebre*DCOS(x) - tht_ebim*DSIN(x))
    q(k) = epsilonz*(qre*DCOS(x) - qim*DSIN(x))
    hs(k) = epsilonz*(hsre*DCOS(x) - hsim*DSIN(x))
    hc(k) = epsilonz*(hcre*DCOS(x) - hcim*DSIN(x))
    WRITE(22,100) u1(k),u2(k),theta1(k),theta2(k), theta_eb(k),q(k) ,hs(k)
  END DO
END DO
GO TO 23
22   CONTINUE

n1=n/k


p= 40000000.d0/l/k
dx=p/n1


DO i=1,k
  
  DO j=1,n1
    x= ((i-1)*n1+j)*dx
    
    IF (x > (40000000.d0/l/2 - p/2) .AND. x < (40000000.d0/l/2 + p/2)) THEN
      
      epsilonz = epsilon
    ELSE
      epsilonz = 0.d0
    END IF
    x= (j-1)*dx*2*pi/p
    k=n1*(i-1)+j
    u1(k) = epsilonz*(u1re*DCOS(x) - u1im*DSIN(x))
    u2(k) = epsilonz*(u2re*DCOS(x) - u2im*DSIN(x))
    theta1(k) = epsilonz*(tht1re*DCOS(x) - tht1im*DSIN(x))
    theta2(k) = epsilonz*(tht2re*DCOS(x) - tht2im*DSIN(x))
    theta_eb(k) = epsilonz*(tht_ebre*DCOS(x) - tht_ebim*DSIN(x))
    q(k) = epsilonz*(qre*DCOS(x) - qim*DSIN(x))
    hs(k) = epsilonz*(hsre*DCOS(x) - hsim*DSIN(x))
    hc(k) = epsilonz*(hcre*DCOS(x) - hcim*DSIN(x))
    WRITE(22,100) u1(k),u2(k),theta1(k),theta2(k), theta_eb(k),q(k),hs(k)
  END DO
END DO
23   CONTINUE
GO TO 40
30   CONTINUE

!      READ*, FILE_NAME
file_name='doutput01'
!      print*,'ENTER DOMAIN SIZE (KM): '
!      READ*,P
p=40000
p= p*1000/l
dx=p/n

OPEN(UNIT=21,FILE=file_name,STATUS='OLD')


DO i=1,n
  
  READ(21,*) u1(i),u2(i),theta1(i),theta2(i),theta_eb(i) ,q(i),hs(i),hc(i)
END DO


40   CONTINUE
101  FORMAT(e25.15,1X,e25.15)
100  FORMAT(e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X  &
    ,e25.15,1X,e25.15,1X,e25.15,1X,e25.15)
RETURN
END SUBROUTINE initial_data

