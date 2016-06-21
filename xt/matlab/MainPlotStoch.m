global xi_s xi_c diurnal

DATADIR = 'data/Ns_100.00-tend_300.00';

ld=1;
RMSP=1; %RMS
AVGP=1; %zonal average
WAVEP=1; % waves
avgHD=1; %high precision average

WAVEPMAVG=0; %waves with avg removed
WAVEFORM=1; %wave format

ebem=11;
xi_s=0.4;
xi_c=0;
tres=3;% hours must be integer.
xgrid=1040; %number of points for the domain
xres=40000/1040;%km
speedres=xres*1000/(3600*6);%speed resolution m/s;
speed=0 ; %m/s
%speed=0;
if (speed==463)
    diurnal=true;
else
    diurnal=false;
end

shift=round(speed/speedres);
snapX=[0,40];
snapT=[400,500];
Xshift=0;

      dels=(16-3*pi*xi_s)/(16+3*pi*xi_s);
      delc=(16-3*pi*xi_c)/(16+3*pi*xi_c);
  
if (ld)
      out=load(fullfile(DATADIR, 'snap_shots'));
end
      
      
      A=load(fullfile(DATADIR, 'rms_energy'));
timet=max(size(A));
timet=(1:1:timet);

%out=out(2000*300+1:end,:);
time=max(size(out))/xgrid;



x=linspace(0,40,xgrid);
pv1=zeros(xgrid,time);
pv2=pv1;pt1=pv1;pt2=pv1;pteb=pv1;pq=pv1;phd=pv1;phc=pv1;phs=pv1;pfc=pv1;
pfd=pv1;pfs=pfd;
    
    
    %moving refrence frame


time=1:1:time*1;
time=time*tres/24;
time=time+timet(end)*tres/24-time(end);
j=ones(size(time'));


RCE=load('SRCE.txt');
fceq=RCE(9);
fdeq=RCE(10);
fseq=RCE(11);




    pv1(:)=circshift(out(:,1),[0,0]);
    pv2(:)=circshift(out(:,2),[0,0]);
    pt1(:)=circshift(out(:,3),[0,0]);
    pt2(:)=circshift(out(:,4),[0,0]) ;
    pteb(:)=circshift(out(:,5),[0,0]);
    pq(:)=circshift(out(:,6),[0,0]) ;
    phs(:)=circshift(out(:,7),[0,0]);
    phc(:)=circshift(out(:,8),[0,0]) ;
    phd(:)=circshift(out(:,9),[0,0]);
    pfc(:)=circshift(out(:,10),[0,0]) ;
    pfd(:)=circshift(out(:,11),[0,0]);
    pfs(:)=circshift(out(:,12),[0,0]);
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
if(avgHD)    
    
outAVG=load('time_aver_out');





    AVGpv1=circshift(outAVG(:,1),[0,0]);
    AVGpv2=circshift(outAVG(:,2),[0,0]);
    AVGpt1=circshift(outAVG(:,3),[0,0]);
    AVGpt2=circshift(outAVG(:,4),[0,0]);
   AVGpteb=circshift(outAVG(:,5),[0,0]);
    AVGpq= circshift(outAVG(:,6),[0,0]);
    AVGphs=circshift(outAVG(:,7),[0,0]);
    AVGphc=circshift(outAVG(:,8),[0,0]);
    AVGphd=circshift(outAVG(:,9),[0,0]);
    AVGpfc=circshift(outAVG(:,10),[0,0]);
    AVGpfd=circshift(outAVG(:,11),[0,0]);
    AVGpfs=circshift(outAVG(:,12),[0,0]);

else

    AVGpv1=mean(pv1');
    AVGpv2=mean(pv2');
    AVGpt1=mean(pt1');
    AVGpt2=mean(pt2');
    AVGpteb=mean(pteb');
    AVGpq=mean(pq');
    AVGphs=mean(phs');
    AVGphc=mean(phc');
    AVGphd=mean(phd');
    AVGpfc=mean(pfc');
    AVGpfd=mean(pfd');
    AVGpfs=mean(pfs');
end
    

    
    
    
meanP=mean(mean(phd+xi_c*phc+xi_s*phs));
meanfc=mean(mean(phc*xi_c))/meanP
meanfs=mean(mean(phs*xi_s))/meanP
meanP=meanP*15.30612244/( 8.1509255/24)
 %   phd=min(phd,.5);
  %  phc=min(phc,.05);




    
     
    
    




if(RMSP)
figure(10)
plot(timet*tres/24,A)
title('RMS of all variables (in nondimentional units)');
legend('$v_1$','$v_2$','$\theta_1$','$\theta_2$','$\theta_{eb}$','q','$h_s$','$h_c$')
xlabel('Days')
 saveas(gcf, 'rms', 'fig')
 print(gcf, '-dpdf', '-r300', 'rms.pdf')

 

end








xindex=[find(x>=snapX(1),1,'first'),find(x<=snapX(2),1,'last')];
tindex=[find(time>=snapT(1),1,'first'),find(time<=snapT(2),1,'last')];

xindex=xindex(1):1:xindex(end);
tindex=tindex(1):1:tindex(end);

x=x(find(x>=snapX(1)));
x=x(find(x<=snapX(2)));
time=time(find(time>=snapT(1)));
time=time(find(time<=snapT(2)));

pv1=pv1(xindex,:);
pv1=pv1(:,tindex);


pv2=pv2(xindex,:);
pv2=pv2(:,tindex);

pt1=pt1(xindex,:);
pt1=pt1(:,tindex);

pt2=pt2(xindex,:);
pt2=pt2(:,tindex);

pteb=pteb(xindex,:);
pteb=pteb(:,tindex);

pq=pq(xindex,:);
pq=pq(:,tindex);


phs=phs(xindex,:);
phs=phs(:,tindex);


phc=phc(xindex,:);
phc=phc(:,tindex);


phd=phd(xindex,:);
phd=phd(:,tindex);

pfc=pfc(xindex,:);
pfc=pfc(:,tindex);


pfd=pfd(xindex,:);
pfd=pfd(:,tindex);

pfs=pfs(xindex,:);
pfs=pfs(:,tindex);





%columnloc=11;
%'spatial location in KM'
%columnloc*xres




%columnplotXT


if(WAVEFORM==1)
    if(AVGP)
    plotAVGHD
    %plotAVGHD2
    end
%plotfieldjournal
%plotfieldjournalsplit

    
    if(WAVEP)
    plotfields3stoch
    end
 
    if(WAVEPMAVG)
    plotfields3MTAVG
    end
else
        if(WAVEP)
    plotfields3stoch2
    end
    if(AVGP)
    plotAVG3stoch3
    end
    if(WAVEPMAVG)
    plotfields3stochMavg2
    end
    
end
                                                                     
