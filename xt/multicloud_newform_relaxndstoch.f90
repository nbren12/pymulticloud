
! Code converted using TO_F90 by Alan Miller
! Date: 2014-04-30  Time: 02:18:49

PROGRAM three_cloud_types_new_formulation
use stochastic
use util
use forcings
use nonlinear_module
use state_mod, only: n, ntrunc, fcls, fdls, fsls, scmt
use cmt_mod
IMPLICIT NONE
INTEGER :: i, j, niter, iter,k,nout, iout,inrg
integer :: GAUGE_ID = 37, SNAPSHOT_ID = 38,ufid=39, tfid=40, PARAMETER_OUTPUT=35
logical :: BINARY_OUTPUT
REAL*8 taumult ,nstochloc
PARAMETER(nstochloc=30)
PARAMETER(taumult=1.d0)

! Toggles
logical, parameter :: CMT_ON = .true.

!     PROGNOSTIC VARIABLES
REAL*8 u1(n),u2(n),theta1(n),theta2(n),  &
    q(n),theta_eb(n),hd(n),hs(n),hc(n),hds(n)
REAL*8 uc(2*ntrunc+1,-1:n+2)
REAL*8 thteb_st(n), asst,lsst






REAL*8 q_tld,alpha_tld,lmd_tld,tempg,dctime,dcint,s1

PARAMETER(alpha_tld=0.1D0,lmd_tld=0.8D0,q_tld=0.9D0)

REAL*8 dt,dx,p, umax,time,dt_max,twave_count2
REAL*8 theta_d,tend,tout,tenergy





! Parameters for precipitation relaxation scheme

REAL*8 fs,fd,fc,deltas
!      PARAMETER(a1=.9,a2=.1,a0=1.7,alpha2=1.d0)



REAL*8 pi,n2,gamma,gammam, cp,omega,r, theta0,g, beta, lcp,EQ
REAL*8 zm, zp, u0, cd
REAL*8 theta_embar

REAL*8 theta_ebs, theta_ebbar

REAL*8 u1_av(n),u2_av(n),tht1_av(n),tht2_av(n),tht_eb_av(n),  &
    q_av(n),hs_av(n),hc_av(n),twave_count,twave_out, pr0_av(n),hd_av(n),  &
    lambda_av(n),lambda,tht_eb_tht_em,pr0


!     on fly averaging

REAL*8 u1a(n),u2a(n),theta1a(n),theta2a(n),theta_eba(n)
REAL*8 qa(n),hsa(n),hca(n),hda(n),fclsa(n),fdlsa(n),fslsa(n)
REAL*8 twave_outa,twave_count2a,tenergya,twa

namelist /data/ tend, tenergy, asst, toggle_nonlinear, stochastic_cmt



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
gamma=1.7D0    ! ratio of moist and dry lapse rates
gammam = 6.0 / km ! moist lapse rate (K /m)
theta0 = 300  ! reference pot temp (K)
omega= 2*pi/(24*hour)! the angular velocity of the earth and
r=6378*km           ! radius of the earth
theta0=300.d0 !K background  reference temperature
g=9.8D0 !m/s gravitational acceleration
beta=2*omega/r      !1/s/m
EQ=40000*km!  Earth's peremeter at the Equator
n2=.0001D0         !Vaissala buoyancy frequency squared


zt=16.00*km !topospheric height
zm = 5000.d0 !middle trop height
zb=500.d0 !meters boundary layer depth
zp = 8.d0 * km ! average height of penetrative clouds

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


! Stochastic Params

CAPE0 = 400 ! J / Kg
MOIST0 = 30 ! K
rstoch = 2 * zp* cp * gammam / theta0 *alpha_bar / c /c

OPEN(35,FILE="OUTPUT",STATUS='new')
WRITE(35,*)'Evaporative time tau_e=', tau_e/3600/24,' days'

PRINT*,'Evaporative time tau_e=', tau_e/3600/24,' days'
WRITE(35,*) 'Type of RCE:'
WRITE(35,*) '-------------'
IF(theta_eb_m_theta_em <= thetam) THEN
  lambdabar=lambdas
  
  WRITE(35,*)'Pure deep convective RCE'
ELSE IF(theta_eb_m_theta_em < thetap) THEN
  WRITE(35,*)'Mixed Deep convective-Congestus RCE'
  lambdabar=a*theta_eb_m_theta_em+b
  
ELSE
  WRITE(35,*)'Pure Congestus RCE--dry troposphere'
  lambdabar=1.d0
END IF

WRITE(35,*)'********lambda_tilde=', lmd_tld


WRITE (*,*) 'QR01'
WRITE (*,*) qr01

WRITE(35,*)'RCE values'
WRITE(35,*)'QR01=', qr01*day, ' K/day Qc=', qcbar*day,' K/day,'  &
    ,'   D=',dbar, 'K m/s', 'QR02=', qr02*day,' K/day'

WRITE(35,*)' theta_eb - theta_em', theta_eb_m_theta_em, 'K'
PRINT*, 'theta_ebs_m_theta_eb', theta_ebs_m_theta_eb, ' K'
print *, 'N Stoch GL' , nstochgl


WRITE(35,*)'-----------------Type of Parametrization'

WRITE(35,*)'a1=',a1,' a2=', a2,' a0=',a0
PRINT*, 'a1=',a1,' a2=', a2,' a0=',a0
WRITE(35,*)'Lambda_bar=', lambdabar
PRINT*,'Lambda_bar=', lambdabar
WRITE(35,*)'mu=', mu
PRINT*,'mu=', mu
WRITE(35,*)'alpha_c=',alpha_c,' tau_c=', tau_c/3600,  &
    ' hours; tau_s=', tau_s/3600, ' hours;tau_conv=', tau_conv/3600,  &
    ' hours;'

PRINT*,'alpha_c=',alpha_c,' tau_c=', tau_c/3600,  &
    ' hours; tau_s=', tau_s/3600, ' hours;tau_conv=', tau_conv/3600,  &
    ' hours;'




!  Reference scales

 ! m isture: ratio of latent heat of vaporization and heat capacity
!times (p/p_0)^{-.28), p=700 hPa, P_0=1000 normalized by alpha_bar unit of temperature.



!    NON DIMENTIONALIZATION OF RCE VALUES AND PARAMETERS USED BEYOND THIS

!    stochastic RCE parameters
! TODO : refactor m0, rstoch, alpha_Bar, l, c into calculate_rce

! FMK13 Taus
tau01 = 1.0d0
tau02 = 3.0d0
tau10 = 1.0d0
tau12 = 1.0d0
tau20 = 3.0d0
tau30 = 5.0d0
tau23 = 3.0d0

! C_omega times from Peters et al. 2013
!tau01 = 1.0d0
!tau02 = 2.2d0
!tau10 = 1.2d0
!tau12 = 1.2d0
!tau20 = 2.4d0
!tau30 = 4.0d0
!tau23 = 0.16d0

if (.false.) then

        OPEN(69,FILE="SRCE.txt",STATUS='old')
        WRITE (*,*) 'rcefile'
        READ(69,*) qbar
        READ(69,*) dbar
        READ(69,*) pbar
        READ(69,*) qr02
        READ(69,*)  tau_e
        READ(69,*)  m0
        READ(69,*) cape0!done
        READ(69,*) capebar
        READ(69,*) fceq!done
        READ(69,*) fdeq!done
        READ(69,*) fseq!done
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
        READ(69,*) moist0!done
endif


tau01=tau01*taumult
tau02=tau02*taumult
tau10=tau10*taumult
tau12=tau12*taumult
tau20=tau20*taumult
tau30=tau30*taumult
tau23= tau23*taumult


call nondimensionalize_params
call init_cmt(n, ntrunc)
call calculate_rce()
! call check_srce()
print *, 'Equilibrium Cloud Fractions:  ', fceq, fdeq, fseq


DO i=1,n
  fcls(i)=fceq;
  fdls(i)=fdeq;
  fsls(i)=fseq;
END DO


WRITE(*,*) 'tau10'
WRITE(*,*) tau10
PRINT*,'m0',m0, 'mu/Qc',mu
WRITE(35,*)'m0',m0*c,'    m/s'
WRITE(PARAMETER_OUTPUT, *) 'Equilibrium Cloud Fractions:  ', fceq, fdeq, fseq
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                  Setup Time Stepping and Output                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

BINARY_OUTPUT = .true.

dt_max =2.d0*minute/t
dt_max =1.d0*minute/t




! tend =  200.0d0 *day/t

niter=tend/dt_max
nout=1

tout=tend/nout

tenergya= dt_max*3.d0
twave_outa=tend-100.d0*day/t

tenergy= (3.d0)*hour/t
! twave_out= 000.d0*day/t



iout=0
inrg=0
OPEN(UNIT=12,FILE='max_amp',STATUS='unknown')
OPEN(UNIT=14,FILE='min_amp',STATUS='unknown')
OPEN(UNIT=16,FILE='rms_energy',STATUS='unknown')
OPEN(UNIT=18,FILE='time_aver_out',STATUS='unknown')
!OPEN(UNIT=SNAPSHOT_ID,FILE='snap_shots',STATUS='unknown')
OPEN(UNIT=11,FILE='time_of_time_aver_out',STATUS='unknown')

write(unit=snapshot_id) n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                          Initialization                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


uc = 0.0d0
CALL initial_data(u1,u2,theta1,theta2,theta_eb,q,hs,hc,n,dx,l,p)
!(u1,u2,theta1,theta2,theta_eb,q,hs,hc,N,DX,L,P)
scmt = 1


WRITE(35,*) 'DX=',dx
WRITE(35,*)'P=',p
WRITE(35,*)'DX*L=',dx*l,' P*L=',p*l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                   SST and Diurnal Cycle Setup                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Prescribed STT with non-zero gradient

!        ASST : strenght of warm pool heating

!       LSST : width of warm pool.

asst=0.5d0;
lsst=EQ/l/8;
!      Strength of the diurnal cycle in K
dctime=0.d0*day/t;
dcint=0.01D0*day/t;
sdc=0.0D0/theta_d;
dtype=2;



! open(8, file='input.nml', status='new')
! write(8,nml=data)
! close(8)
open(8, file='input.nml', status='old')
read(8,nml=data)

tend = tend * day / t
tenergy = tenergy * hour/t

print *, 'data output frequency for waves:', tenergy*t/hour,' hours'
print *, 'total run time', tend*t/day,'days'
print *, 'asst=', asst
print *, 'toggle_nonlinear=', toggle_nonlinear
print *, 'stochastic_cmt=', stochastic_cmt
close(8)
!        SEPT in BD (non dim units)

OPEN(UNIT=36,FILE='sst',STATUS='unknown')

DO i=1,n
  IF(i == n)PRINT*,i*dx,EQ/l/2,2*lsst,i*dx-EQ/l/2-2*lsst
  IF(DABS(i*dx-EQ/l/2) <= 2*lsst) THEN
    thteb_st(i)=asst*theta_ebs_m_theta_eb *COS(pi*(i*dx-EQ/l/2)/2/lsst)
  ELSE
    thteb_st(i)=-asst*theta_ebs_m_theta_eb;
  END IF
  WRITE(36,*) thteb_st(i)
END DO



CLOSE(36)
!       STOP
WRITE(35,*)'DBAR "non_dim"', m0*(1.d0 - mu*qr02)*lambdabar  &
    *theta_eb_m_theta_em, 'DBAR=',dbar
WRITE(35,*)'UD:', ud
WRITE(35,*)'A=', a, 'B=', b, 'mu =',mu



!     output initial data

time =0.d0
CALL output(u1,u2,theta1,theta2,theta_eb,q,hs,hc,hd,n,iout)
CALL out_amp(u1,u2,theta1,theta2,theta_eb,q,hs,hc,n,time)
WRITE(35,*) 'BACK to MAIN'


twave_count= 0.d0
DO i=1,n
  
  u1_av(i) = 0.d0
  u2_av(i) = 0.d0
  tht1_av(i) = 0.d0
  tht2_av(i) = 0.d0
  tht_eb_av(i) = 0.d0
  q_av(i) = 0.d0
  hs_av(i) = 0.d0
  hc_av(i) = 0.d0
  pr0_av(i)=0.d0
  hd_av(i)=0.d0
  lambda_av(i)=0.d0
  hds(i)=0.d0
  hs(i)=0.d0
  u1a(i)=0.d0
  u2a(i)=0.d0
  theta1a(i)=0.d0
  theta2a(i)=0.d0
  theta_eba(i)=0.d0
  qa(i)=0.d0
  hsa(i)=0.d0
  hca(i)=0.d0
  hda(i)=0.d0
  fclsa(i)=0.d0
  fdlsa(i)=0.d0
  fslsa(i)=0.d0
  
END DO





twave_count2a=0.d0
twave_count2=0.d0
iter=0


open(unit=GAUGE_ID, file='gauge')

if ( BINARY_OUTPUT) then
   OPEN(UNIT=SNAPSHOT_ID, FILE = 'snap_shots',&
        FORM = 'unformatted', access='stream', status="NEW")

   open(unit=18, file='real.bin', form='unformatted', status='NEW', access='stream')
   write (18) n
   write (18) ntrunc
   write (18) dx*L
   write (18) T
   write (18) L
   write (18) alpha_bar
   close(18)
else
   OPEN(UNIT=SNAPSHOT_ID,FILE='snap_shots',STATUS='unknown')
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         Begin Iterations                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


WRITE(35,*)'BEGIN ITERATIONS:'
WRITE(35,*)'********************************************'
10   CONTINUE
iter=iter+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               Select dt according to CFL condition                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

umax=1.d0
DO i=1,n
  umax = DMAX1(umax, (DABS(q(i)+q_tld)+DABS( alpha_tld*q(i) +  &
      lmd_tld*q_tld)+ DABS(u1(i)+alpha_tld*u2(i)))/5.d0)
END DO




dt = DMIN1(0.9*dx/umax,dt_max)/2d0


! Abort conditions
if (dt < 1d0/t) then
   print *, 'Timestep too small'
   print *, 'Aborting run!!!!'
   stop 9
end if
   do i=1,n
      if (uc(1,i) /= uc(1,i)) then
         print *, 'NA in velocity field'
         print *, 'Aborting run!!!!'

         stop 9
      end if
   end do

IF(iter == 1)THEN
  WRITE(35,*)' Initial Time step =', 2*dt*t/minute, 'minutes  iter=',     iter
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                       Step Forward in Time                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

i = 50
WRITE(GAUGE_ID, 103) time * t / day,c* u1(i),c*u2(i), alpha_bar * theta1(i),alpha_bar * theta2(i),  &
    alpha_bar * theta_eb(i), q(i), &
    alpha_bar/t *day* DMAX1(0.00000000000001D0,hs(i)),  &
    alpha_bar/t *day*DMAX1(0.00000000000001D0,hc(i)), &
    alpha_bar/t *day* DMAX1(0.00000000000001D0,hd(i)),  &
    fcls(i),fdls(i),fsls(i)

! unpack data from u vector
! adjust temperature from FMK13 convection : m  theta_m -> theta_m
DO i=1,n
      uc(1,i) = u1(i)
      uc(2,i) = u2(i)
      uc(1+ntrunc,i) = theta1(i)*1
      uc(2+ntrunc,i) = theta2(i)*2
      uc(2*ntrunc+1,i) = q(i)
END DO

CALL central_scheme(uc,dx,dt,n,q_tld,alpha_tld,lmd_tld, ntrunc)

call vertical_advection_driver(uc, dx, 2d0*dt, n, ntrunc)

if (.not. cmt_on) then
   ! damping (comment beacuse it is inside cmt)
   do i=1,n
      do j = 1,ntrunc
         uc(j,i) = dexp(-ud*2d0*dt) * uc(j,i)
      end do
   end do

end if

! split-in-time vertical advection

! copy output to flux array for stochastic multicloud step
! readdjust temperature to FMK13 convection :  theta_m -> m theta_m
DO i=1,n
  u1(i) = uc(1,i)
  u2(i) = uc(2,i)
  theta1(i) = uc(ntrunc+1,i)/1
  theta2(i) = uc(ntrunc+2,i)/2
  q(i) = uc(2*ntrunc+1,i)
END DO

IF (time > dctime) THEN
  s1=DMIN1(1.d0,(time-dctime)/dcint)
END IF




tempg=dt
CALL updatehcds(fcls,fdls,fsls,u1, u2, theta1,theta2,theta_eb,q,hds,hc,hd  &
    ,n, dx, 2.d0*dt*t/(hour) )
dt=tempg

CALL range_kuttas(u1,u2,theta1,theta2,theta_eb,q,hs,hc,hd,n  &
    ,2.d0*dt,thteb_st,time,hds,s1)
dt=tempg



! unpack data from u vector
! adjust temperature from FMK13 convection : m  theta_m -> theta_m
DO i=1,n
   uc(1,i) = u1(i)
   uc(2,i) = u2(i)
   uc(1+ntrunc,i) = theta1(i)*1
   uc(2+ntrunc,i) = theta2(i)*2
   uc(2*ntrunc+1,i) = q(i)
END DO

if (cmt_on) call updatecmt(uc(1:ntrunc,1:n), scmt, hd, hc, hs, 2d0*dt)
CALL central_scheme(uc,dx,dt,n,q_tld,alpha_tld,lmd_tld, ntrunc)

! copy output to flux array for stochastic multicloud step
! readdjust temperature to FMK13 convection :  theta_m -> m theta_m
DO i=1,n
   u1(i) = uc(1,i)
   u2(i) = uc(2,i)
   theta1(i) = uc(ntrunc+1,i)/1
   theta2(i) = uc(ntrunc+2,i)/2
   q(i) = uc(2*ntrunc+1,i)
END DO

time= time+2*dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                      Output Amplitude Data                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(time >= (inrg+1)*tenergy) THEN
  CALL out_amp(u1,u2,theta1,theta2,theta_eb,q,hs,hc,n,time)
  inrg=inrg+1
!         write(*,*) inrg
  
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         Snapshot Output                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF(time >= (twave_count2+1)*tenergy) THEN
  
  twave_count2=twave_count2+1
  
  if (BINARY_OUTPUT) then

     print *,  'outputing snapshot at time', time * T / day
     ! WRITE(unit=SNAPSHOT_ID) time *T / day, c*u1, c*u2,&
     !      alpha_bar *theta1, alpha_bar * theta2,  &
     !      alpha_bar*theta_eb, &
     !      alpha_bar * q,&
     !      alpha_bar/ (T / day) *DMAX1(0D0,hs),  &
     !      alpha_bar/ (T / day) *DMAX1(0D0,hc), &
     !      alpha_bar/ (T / day) *DMAX1(0D0,hd),  &
     !      fcls,fdls,fsls

     open(unit=ufid, file='real.bin', position='append', status='unknown', form='unformatted', access='stream')
     write(ufid) time * T/ day
     write(ufid) uc(1:ntrunc,1:n) * c
     write(ufid) uc(ntrunc+1:2*ntrunc,1:n) * alpha_bar
     write(ufid) uc(2*ntrunc+1,1:n) * alpha_bar
     write(ufid) theta_eb * alpha_bar
     write(ufid) hc * alpha_bar/ (T / day)
     write(ufid) hd * alpha_bar/ (T / day)
     write(ufid) hs * alpha_bar/ (T / day)
     write(ufid) fcls
     write(ufid) fdls
     write(ufid) fsls
     write(ufid) scmt
     close(ufid)

  else
      DO i=1,n
    

          WRITE(SNAPSHOT_ID,102) u1(i),u2(i), theta1(i),theta2(i),  &
              theta_eb(i), q(i),DMAX1(0.00000000000001D0,hs(i)),  &
              DMAX1(0.00000000000001D0,hc(i)), DMAX1(0.00000000000001D0,hd(i)),  &
              fcls(i),fdls(i),fsls(i)
    
      END DO
  !               write(*,*) twave_count2
  endif
  
  
  
END IF




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                     Update Time Average data                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF((time+twave_outa-tend) >= (twave_count2a+1)*tenergya) THEN
  
  twave_count2a=twave_count2a+1
  twa=twave_count2a
  
  DO i=1,n
    
    u1a(i)=u1a(i)+u1(i)
    u2a(i)=u2a(i)+u2(i)
    theta1a(i)=theta1a(i)+theta1(i)
    theta2a(i)=theta2a(i)+theta2(i)
    theta_eba(i)=theta_eba(i)+theta_eb(i)
    qa(i)=qa(i)+q(i)
    hsa(i)=hsa(i)+hs(i)
    hca(i)=hca(i)+hc(i)
    hda(i)=hda(i)+hd(i)
    fclsa(i)=fclsa(i)+fcls(i)
    fdlsa(i)=fdlsa(i)+fdls(i)
    fslsa(i)=fslsa(i)+fsls(i)
  END DO
  
  
  
END IF



! IF(time >= (iout+1)*tout) THEN
  
!   WRITE(35,*)'Time step=', 2*dt*t/minute,' minutes',' iter=',iter
!   iout=iout+1
!   WRITE(35,*)'iout=',iout,' time=',time*t/hour,' hours'
!   CALL output(u1,u2,theta1,theta2,theta_eb,q,hs,hc,hd,n,iout)
  
! END IF
IF(time < tend) GO TO 10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                    End Iterations...Finish up                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




WRITE(35,*)'*************************************************'
WRITE(35,*)'Program ended.Output Final results:',iout,' time=',  &
    time*t/hour,' hours'

CALL output(u1,u2,theta1,theta2,theta_eb,q,hs,hc,hd,n,iout)
WRITE(35,*)'***********************END**************************'


!     save avg

DO i=1,n
  
  WRITE(18,102) u1a(i)/twa,u2a(i)/twa, theta1a(i)/twa,theta2a(i)/twa,  &
      theta_eba(i)/twa, qa(i)/twa,hsa(i)/twa,hca(i)/twa,hda(i)/twa,  &
      fclsa(i)/twa,fdlsa(i)/twa,fslsa(i)/twa
END DO



CLOSE(12)
CLOSE(14)
CLOSE(16)
CLOSE(18)
CLOSE(11)
close(GAUGE_ID)
close(SNAPSHOT_ID)

100  FORMAT(e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X,e25.15)
101  FORMAT(e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X,e25.15,  &
    1X,e25.15,1X,e25.15)



102  FORMAT(e25.15,1X,e25.15,1X, e25.15,1X,e25.15,1X,  &
    e25.15,1X e25.15,1X,e25.15,1X,  &
    e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X)
103  FORMAT(e25.15,e25.15,1X,e25.15,1X, e25.15,1X,e25.15,1X,  &
    e25.15,1X e25.15,1X,e25.15,1X,  &
    e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X)

CONTAINS


subroutine nondimensionalize_params()



    t=l/c            !~ 8.33 hours, time scale
    theta_d=alpha_bar !~15 degree K; temperature scale

    WRITE(35,*)'L=',l,' meters,  T=',t,'sec,  C=',c,  &
        ' met./sec; theta_d=alpha_bar=',alpha_bar, 'Kelvin'


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
    WRITE (*,*) 'zb'
    WRITE (*,*) zb

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

    print *, 'Should be close to zero'
    print *,fvec
    print *
end subroutine

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

end subroutine

END PROGRAM three_cloud_types_new_formulation


!     MAIN Program Ends here




SUBROUTINE out_amp(u1,u2,theta1,theta2,theta_eb,q,hs,hc,n,time)
IMPLICIT NONE

REAL*8, INTENT(IN)                       :: u1(n)
REAL*8, INTENT(IN)                       :: u2(n)
REAL*8, INTENT(IN)                       :: theta1(n)
REAL*8, INTENT(IN)                       :: theta2(n)
REAL*8, INTENT(IN)                       :: theta_eb(n)
REAL*8, INTENT(IN)                       :: q(n)
REAL*8, INTENT(IN)                       :: hs(n)
REAL*8, INTENT(IN)                       :: hc(n)
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN OUT)                   :: time
INTEGER :: i


REAL*8 mxu1,mxu2,mxtht1,mxtht2,mxtht_eb,mxq,mxhs,mxhc
REAL*8 mnu1,mnu2,mntht1,mntht2,mntht_eb,mnq,mnhs,mnhc
REAL*8 rmsu1,rmsu2,rmstht1,rmstht2,rmstht_eb,rmsq,rmshs,rmshc

!         mxu1=0.d0
!         mnu1=0.d0
!         mxu2=0.d0
!         mnu2=0.d0
!         mxtht1=0.d0
!         mxtht2=0.d0
!         mntht1=0.d0
!         mntht2=0.d0
!         mxtht_eb=0.d0
!         mntht_eb=0.d0
!         mxq=0.d0
!         mnq=0.d0
!         mxhs=0.d0

!         mnhs=0.d0
!         mxhc=0.d0

!         mnhc=0.d0

rmsu1=0.d0
rmsu2=0.d0
rmstht1=0.d0
rmstht2=0.d0
rmshs=0.d0
rmshc=0.d0

rmstht_eb=0.d0
rmsq=0.d0

DO i=1,n
!            mxu1=dmax1(mxu1,u1(i))
!            mxu2=dmax1(mxu2,u2(i))
!            mxtht1=dmax1(mxtht1,theta1(i))
!            mxtht2=dmax1(mxtht2,theta2(i))
!            mxtht_eb=dmax1(mxtht_eb,theta_eb(i))
!            mxq=dmax1(mxq,q(i))
!            mxhs=dmax1(mxhs,hs(i))
!            mxhc=dmax1(mxhc,hc(i))
  
!            mnu1=dmin1(mnu1,u1(i))
!            mnu2=dmin1(mnu2,u2(i))
!            mntht1=dmin1(mntht1,theta1(i))
!            mntht2=dmin1(mntht2,theta2(i))
!            mntht_eb=dmin1(mntht_eb,theta_eb(i))
!            mnq=dmin1(mnq,q(i))
!            mnhs=dmin1(mnhs,hs(i))
!            mnhc=dmin1(mnhc,hc(i))
  
  
  rmsu1=rmsu1+u1(i)**2/n
  rmsu2=rmsu2+u2(i)**2/n
  rmstht1=rmstht1+theta1(i)**2/n
  rmstht2=rmstht2+theta2(i)**2/n
  rmstht_eb=rmstht_eb+theta_eb(i)**2/n
  rmsq=rmsq+q(i)**2/n
  rmshs=rmshs+hs(i)**2/n
  rmshc=rmshc+hc(i)**2/n
  
  
END DO



!         write(12,101)  mxu1,mxu2,mxtht1,mxtht2,mxtht_eb,mxq,mxhs,mxhc
!         write(14,101)  mnu1,mnu2,mntht1,mntht2,mntht_eb,mnq,mnhs,mnhc
WRITE(16,101)   DSQRT(rmsu1),DSQRT(rmsu2),DSQRT(rmstht1),  &
    DSQRT(rmstht2),DSQRT(rmstht_eb),DSQRT(rmsq), DSQRT(rmshs),DSQRT(rmshc)

101     FORMAT(e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X,e25.15  &
    ,1X,e25.15,1X,e25.15)
RETURN
END SUBROUTINE out_amp


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
READ*, idata
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
WRITE(35,*)'Wavelength = ', p*l/1000.d0, ' km'
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
WRITE(35,*)'Sine wave packet'
kn=k
n1=n/kn


p= 40000000.d0/l/kn
dx=p/n1

WRITE(35,*)' Small Wavelength = ', p*l/1000.d0, ' km'
WRITE(35,*)' Large Wavelength = ', k*p*l/1000.d0, ' km'

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
WRITE(35,*)'localized single wave'

n1=n/k


p= 40000000.d0/l/k
dx=p/n1

WRITE(35,*)' Small Wavelength = ', p*l/1000.d0, ' km'
WRITE(35,*)' Large Wavelength = ', k*p*l/1000.d0, ' km'

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
!      WRITE(35,*)'ENTER NAME OF (PUTPUT) DATA FILE: '
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


SUBROUTINE output(u1,u2,theta1,theta2,theta_eb,q,hs, hc,hd,n,nout)

IMPLICIT NONE
REAL*8, INTENT(IN OUT)                   :: u1(n)
REAL*8, INTENT(IN OUT)                   :: u2(n)
REAL*8, INTENT(IN OUT)                   :: theta1(n)
REAL*8, INTENT(IN OUT)                   :: theta2(n)
REAL*8, INTENT(IN OUT)                   :: theta_eb(n)
REAL*8, INTENT(IN OUT)                   :: q(n)
REAL*8, INTENT(IN OUT)                   :: hs(n)
REAL*8, INTENT(IN OUT)                   :: hc(n)
REAL*8, INTENT(IN OUT)                   :: hd(n)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: nout
INTEGER :: i
REAL*8  dx,pi,l, x
CHARACTER (LEN=2) :: cout 
CHARACTER (LEN=8) filename


WRITE(35,*)'*************output #',nout

OPEN(13,FILE='conv',STATUS='unknown')
WRITE(13,101)nout
CLOSE(13)
OPEN(13,FILE='conv',STATUS='unknown')

READ(13,*) cout
CLOSE(13)
!      write(*,*)'cout=',cout(1:2)

filename='output'//cout
!      filename(1:6)='output'
!      filename(7:8)=cout(1:2)
WRITE(35,*)'filename =', filename

OPEN(UNIT=27,FILE=filename,STATUS='unknown')


DO i=1,n
  WRITE(27,100) u1(i),u2(i),theta1(i),theta2(i),theta_eb(i)  &
      ,q(i),hs(i),hc(i),hd(i)
END DO
WRITE(35,*)'output ended'
CLOSE(27)
101   FORMAT('''',i2.2,'''')
100   FORMAT(e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X,e25.15,1X,e25.15,  &
    1X,e25.15,1X,e25.15,1X,e25.15)
RETURN
END SUBROUTINE output
