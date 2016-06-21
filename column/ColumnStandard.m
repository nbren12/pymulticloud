%Stochastic multicloud model with congestus detrainment moistening and lag
%type stratiform heating closure.  
clear all

NDAYS=100;  %legth of simulation run last NDAYS/2 days are plotted





dcyclestrengh=0;
dcycle=dcyclestrengh;
stoch=true;








global Tau Rc xi_s xi_c alpha_s alpha_c
global tau01
global tau02
global tau10
global tau12
global tau20
global tau30
global tau23
global r23value
global timeset
Cd0=10^(-3)
xi_s=0.4;
N2=0.0001; %s^-2, Vaisala frequency, squared,
Ep=0.9; % prec ipip. effeciency
  hr=3600;
day=24*hr; 
km=1000;
QR0=1./day; % K/day , radiative cooling
  theta0=300; %K, reference temperatureå
    g=9.8;%m/s^2, gravitational acceleration
      gamma=1.7;% ratio of moist and dry lapse rates
Gammam=6.0e-3;% K/m,  moist lapse rate
  cp=1.e3; % J/kg/K, specific heat  at constant pressure 
    H=8000; % m, average height of penetrative clouds
      Hm=5000; %m average middle troposphere level
    hb=500; % m, depth of boundary layer   
%      T=8*hr; % reference time scale
ZT=16000;% height of the troposphere
c=sqrt(N2)*ZT/pi;% ~50 m/sec; gravity wave speed; velocity time scale
beta=2.2804e-11; %1/sec 1/meter Grdient of Coriolis force at equator
Le= sqrt(c/beta); %~ 1500*km; % Rossby deformation radius; length time scale
T=Le/c; % ~ 8.33 hours, Time scale    
alpha_bar=Hm*N2*theta0/g;





%*******************************************************************************************
%key parameters for congestus detrainment moistening

tebsmteb = 10/alpha_bar;% theta*_eb-theta_eb
tebmtembar = 11/alpha_bar;% theta_eb-theta_em
tebmtelbar=0/alpha_bar;% theta_eb-theta_el
deltac = 0;
deltacl= 1;
xi_c= 16/3/pi*(1-deltacl)/(1+deltacl);

%*******************************************************************************************
Cd0=10^(-3)

R=2*H*cp*Gammam/theta0;
R = R*alpha_bar/c^2; %Equivalent to a convective time scale of 70 seconds. Too quick, therefore not used here.
%R = 214113*10^(-4)*(c^2)
%use   tau_conv = 2 hours  % as  in determinestic MC models
tau_conv = 2*3600/T;%
%tau_e=8*hr/T;
%tau_eb= (Ep*hb/Hm);tau_em=(Ep*H/Hm);
tau_R=50*day/T; % Newtonian cooling relaxation time
H = H/Le; hb=hb/Le;ZT = ZT/Le;Hm=Hm/Le;
QR01=QR0*T/alpha_bar;
tau_s = 3*3600/T;%
delay=true;
%%%% Parameters
mu=0.25;
gamma2=0.1;
gamma2p=2;
alpha2=0.1;
alphac=0.1;% alphac=0.05 so that at equilibrium with fdeq=0.02 and fceq=0.009, we have Hdbar/Hcbar=44.44, as in multicloud models.
alpha_c=alphac;

a1=0.5;
a2=0.5;
a0=2;

%*******************************************************************************************
% key parameters 
QCAPE=0; %coefficient of q in CAPE parameterization (a_2^c in the report)
LCAPE=true;% use low level CAPE (instead of CAPE) in computation of congestus creation rate
tauscale=1;  %multiply all the transition time scales by tauscale
alpha_s=.25;%as in deterministic multicloud model, used for the lag type stratiform closure.
%*******************************************************************************************


CAPE0 = 400;% J/kg
CAPE0 = CAPE0/c^2;
MOIST0=30/alpha_bar;
n = 30;
DT0 = 5*60/T;%5 minutes
DT=DT0/12;





%set transition timescales******************************************
Tau=zeros(4,4);
timeset=4;
r23value=1;
switch(timeset)
    case(1)
        %time units, hour
    
    tau01=1*tauscale;%
    tau02=3*tauscale;%
    tau10=1*tauscale;
    tau12=1*tauscale;
    tau20=3*tauscale;
    tau30=5*tauscale; 
    
    tau23=3*tauscale;
     
    
    case(2)
         %time units, hour
    
    tau01 =1;%
    tau02=3;%
    tau10=5;
    tau12=1;
    tau20=5;
    tau30=5;
    
    tau23=1.5;
     
   
    case(3)
    %time units, hour
    
    
    tau01 =3;%
    tau02= 5;%
    tau10=2;
    tau12=2;
    tau20=5;
    tau23=.5;
    tau30=3;

    case(4)

    tau01=1.0;%
    tau02= 2.2;%
    tau10=1.2;
    tau12=1.2;
    tau20=2.4;
    tau23=.16;
    tau30=4.0;

end

Tau(1,2)=tau01*hr;
Tau(2,1)=tau10*hr;
Tau(2,3)=tau12*hr;
Tau(1,3)=tau02*hr;
Tau(3,4)=tau23*hr;
Tau(3,1)=tau20*hr;
Tau(4,1)=tau30*hr;




Tau=Tau/T;
Rc=zeros(4,4);



%


% Initial equilibrium values for filling
% fraction of stochastic cloud cover

% Radiative convective equilibrium solution: RCE

tt =tebmtembar/MOIST0;

tmp = @(x)RCEr7(x,QR01,tt,Hm,CAPE0);



X=fsolve(tmp,[.1 .1 .1 1]);
X=real(X);

CAPEbar=X(4);
setRc(CAPEbar/CAPE0,tt)
%[pi_eq,r01,r02,r10,r20,r23,r12,r30] =  equilibriumdistribution(CAPEbar/CAPE0,tt);

fceq = X(1);
fdeq = X(2);
fseq = X(3);

'approximatios'
% fd=1/(1+ xi_s*alpha_s*Tau(4,1)/Tau(3,4))
% fs=1-fd
% 
% xis=(1-fd)/(fd*alpha_s*Tau(4,1)/Tau(3,4))

%fc=xi_c*alpha_c*fceq/( fdeq+ xi_s*alpha_s*fseq + xi_c*alpha_c*fceq)
%fd=fdeq/( fdeq+ xi_s*alpha_s*fseq + xi_c*alpha_c*fceq)
%fs=xi_s*alpha_s*fseq /( fdeq+ xi_s*alpha_s*fseq + xi_c*alpha_c*fceq)




Hdbar = fdeq*sqrt(CAPEbar)/Hm;
Hsbar = fseq*alpha_s*sqrt(CAPEbar)/Hm;
Hcbar = fceq*alpha_c*sqrt(CAPEbar)/Hm;



QR02 = -Hsbar + Hcbar;
tspi = 2*sqrt(2)/pi;
Pbar = Hdbar+Hsbar*xi_s+Hcbar*xi_c;
m0 = Pbar*2*sqrt(2)*ZT/pi;
m0=m0/((1+mu*(Hsbar-Hcbar)/QR01)*tebmtembar +deltac*sqrt(2)*Hcbar*tebmtelbar/(QR01*pi) );
Dbar=m0*(1+mu*(Hsbar-Hcbar)/QR01);
Ecbar=m0*deltac*sqrt(2)*Hcbar*tebmtelbar/(QR01*pi);

taue=hb*tebsmteb/(Ecbar*tebmtelbar+Dbar*tebmtembar);


if (m0 < 0) 
    m0 
    pause
end


Qbar=sqrt(CAPEbar)/Hm
%% 

fdtest=fdeq;fctest=fceq;fstest=fseq;



% stability_check

% 'Check if stable then hit a key to continue or ^C to stop'
% 
% pause
%%
fc=fceq;fd=fdeq;fs=fseq;

z = ones(n,n);

rand(0)
state = 0*floor(4*rand(n,n));

tht1=rand(1,1)*1.0e-5;tht2=rand(1,1)*1e-5;thteb=rand(1,1)*1e-5;q=rand(1,1)*1e-5;
tht1t=0;tht2t=0;thtebt=0;qt=0;

CAPE = max(CAPEbar +R*(thteb - gamma*(tht1+gamma2*tht2)),0);

 
CAPEc = max(CAPEbar +R*(thteb - gamma*(tht1 + gamma2p*tht2) ),0);

% stoch=input('Enter 1 for stochastic MC 0 otherwise (birth death or deterministic')


tebmtem = tebmtembar + thteb - q - tspi*(tht1+alpha2*tht2);
tebmtel=  tebmtelbar + thteb - 2*q - tspi*(tht1+2*tht2);



Hdp = max( fdeq*sqrt(CAPEbar)/Hm + (1/tau_conv)*(a1*thteb+a2*q - a0*(tht1+gamma2*tht2)),0);
Hd  = fd/fdeq*(Hdp);
%Hd = fd*sqrt(CAPE)/Hm;
if(delay==false)
Hsp = max( fseq*sqrt(CAPEbar)/Hm + (1/tau_conv)*(a1*thteb+a2*q - a0*(tht1+gamma2*tht2)),0);

Hs  = (fs/fseq)*(Hsp)/4;
else
    Hs=Hsbar;
end

%Hs =  fs  *sqrt(CAPE)/Hm/4;
Hc =  fc  *sqrt(CAPEc)/Hm*alphac;
%Hc =  fc  *sqrt(CAPE)/Hm*alphac
D = m0*max((1+mu*(Hs-Hc)/QR01),0)*tebmtem;
Ec=tebmtel*max(m0*deltac*sqrt(2)/(QR01*pi),0)*Hc;

 CAPEs=CAPE;
 Ds=D;

FTHT1m1 = Hd + xi_s *Hs +xi_c *Hc - QR01 - tht1/tau_R;
FTHT2m1 = -Hs + Hc - QR02 - tht2/tau_R;
FTHEBm1 = (tebsmteb - thteb)/taue - (D+Ec)/hb;
FQm1 = - tspi*(Hd + xi_s *Hs +xi_c *Hc)+ (D+Ec)/ZT ;
%pause

FTHT1m2 = Hd + xi_s *Hs +xi_c *Hc - QR01 - tht1/tau_R;
FTHT2m2 = -Hs + Hc - QR02 - tht2/tau_R;
FTHEBm2 = (tebsmteb - thteb)/taue - (D+Ec)/hb;
FQm2 = - tspi*(Hd + xi_s *Hs +xi_c *Hc) + (D+Ec)/ZT ;
if(delay)
FHs1=(alpha_s*fs*Hd/fdeq-Hs)/tau_s;
FHs2=(alpha_s*fs*Hd/fdeq-Hs)/tau_s;
end

%%

if(delay)
sol=[tht1,tht2,thteb,q,fceq,fdeq,fseq,Hs];
else
sol=[tht1,tht2,thteb,q,fceq,fdeq,fseq];

end


tend = NDAYS*24*3600/T;

Niter = round(tend/DT0);


DT0=DT0*T/3600;

Hd=Hdbar;
Hs=Hsbar;
Hc=Hcbar;

p23=[];p30=[];

saveSRCE
%tebsmteb=tebsmteb+5/alpha_bar

for I=1:Niter
I/288%dcycle
tebd=dcyclestrengh*tebsmteb*sin(((I/288)-.25)*2*pi);%dcycle

CAPE = max(CAPEbar +R*(thteb+QCAPE*q - gamma*(tht1+gamma2*tht2)),0);


CAPEc = max(CAPEbar +R*(thteb - gamma*(tht1 + gamma2p*tht2) ),0);


tebmtem = tebmtembar + thteb - q - tspi*(tht1+alpha2*tht2);
tebmtel=  tebmtelbar + thteb - 2*q - tspi*(tht1+2*tht2);





DRY = max(tebmtem,0)/MOIST0;

if(stoch)

  %[fc,fd,fs]=KMMb2(fc,fd,fs,CAPE/CAPE0,DRY,DT, n);
if(LCAPE)
  [fc,fd,fs]=birthdeathexact2fort(fc,fd,fs,CAPE/CAPE0,CAPEc/CAPE0,DRY,DT0, n);
else
  [fc,fd,fs]=birthdeathexact2(fc,fd,fs,CAPE/CAPE0,CAPE/CAPE0,DRY,DT0, n);
end

end

if(~stoch)
    
  [fc,fd,fs]=mf(fc,fd,fs,CAPE/CAPE0,CAPEc/CAPE0,DRY,DT0, n,T);
    
end



Hdp = max( fdeq*sqrt(CAPEbar)/Hm + (1/tau_conv)*(a1*thteb+a2*q - a0*(tht1+gamma2*tht2)),0);
Hd  = (fd/fdeq)*(Hdp);

if(delay==false)
Hsp = max( fseq*sqrt(CAPEbar)/Hm + (1/tau_conv)*(a1*thteb+a2*q - a0*(tht1+gamma2*tht2)),0);
Hs  = (fs/fseq)*(Hsp)/4;
end
%Hs =  fs  *sqrt(CAPE)/Hm/4;
%Hc =  fc  *sqrt(CAPE)/Hm*alphac;



%Hsp = max( fseq*sqrt(CAPEbar)/Hm + (1/tau_conv)*(a1*thteb+a2*q - a0*(tht1+gamma2*tht2)),0);
%Hs  = (fs/fseq)*(Hsp)/4;



Hc =  alphac*fc*sqrt(CAPEc)/Hm;
D = m0*max((1+mu*(Hs-Hc)/QR01),0)*tebmtem;
Ec=tebmtel*max(m0*deltac*sqrt(2)/(QR01*pi),0);

CAPEs=[CAPEs;CAPE];Ds=[Ds;D];

FTHT1 = Hd + xi_s *Hs +xi_c *Hc - QR01 - tht1/tau_R;
FTHT2 = -Hs + Hc - QR02 - tht2/tau_R;
FTHEB = (tebsmteb+tebd - thteb)/taue - (D+Ec)/hb;%dcycle
FQ = - tspi*(Hd + xi_s *Hs +xi_c *Hc) + (D+Ec)/ZT ;
FHs=(alpha_s*fs*Hd/fdeq-Hs)/tau_s;
%pause

tht1 = tht1 + DT*(23*FTHT1-16*FTHT1m1+5*FTHT1m2);

tht2 = tht2 +  DT*(23*FTHT2-16*FTHT2m1+5*FTHT2m2);
thteb = thteb + DT*( 23*FTHEB-16*FTHEBm1+5*FTHEBm2) ;
q = q + DT*(23*FQ-16*FQm1+5*FQm2);
if(delay)
Hs = Hs + DT*(23*FHs-16*FHs1+5*FHs2);
end
if(delay)
sol=[sol; tht1,tht2,thteb,q,fc,fd,fs,Hs];
else
    sol=[sol; tht1,tht2,thteb,q,fc,fd,fs];

end
FTHT1m2 = FTHT1m1;
FTHT2m2 = FTHT2m1;
FTHEBm2 = FTHEBm1;
FQm2 = FQm1;
if(delay)
FHs2 = FHs1;
end
FTHT1m1 = FTHT1;
FTHT2m1 = FTHT2;
FTHEBm1 = FTHEB;
FQm1 = FQ;
if(delay)
FHs1 = FHs;
end

%tht1t=tht1;tht2t=tht2;thtebt=thteb;qt=q;
%tht1=tht1p;tht2=tht2p;thteb=thtebp; q=qp;

%break

end

if(stoch)
        WLINE=2;
else
    WLINE=6;
end
columnplotStoch2
if(dcyclestrengh>0)
   dplot
end








