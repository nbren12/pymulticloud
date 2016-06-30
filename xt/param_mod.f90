module param

  real(8) :: hour, minute, day, km, tempouttime
  parameter(hour=3600.d0, minute=60.0d0, day=86400.0d0, km=1000.d0)

  ! basic nondimensional scales
  real(8) :: t, c, l, alpha_bar


  real(8) nstochgl

  real(8) a,b,a0,a1,a2,a0p,a1p,a2p,alpha2,alpha3,  &
       lambdas,alpha_s,alpha_c,xis,xic,mu,ud,thetap,thetam,zb,zt
  real(8) tau_conv,tau_e, tau_s,tau_c,tau_r, tau_d, sdc,dtype

  real(8) tau01,tau02,tau10,tau12,tau20,tau30,tau23, r23value, times
  real(8) capebar, cape0,dry,pbar, qbar,fceq,fdeq,fseq, rstoch

  real(8) theta_eb_elbar,deltac,deltac1

  real(8) qr01,qr02,lambdabar,qcbar,dbar,  &
       theta_eb_m_theta_em,m0,theta_ebs_m_theta_eb,moist0,alpha4

end module param
