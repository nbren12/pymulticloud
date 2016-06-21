sol=sol(round(max(size(sol)/2)):end,:);
Niter=max(size(sol))-1;

tht1=sol(:,1);
tht2=sol(:,2);
thteb=sol(:,3);
q=sol(:,4);
fc=sol(:,5);
fd=sol(:,6);
if(delay)
fs=sol(:,7);
else
end
CAPE = max(CAPEbar +R*(thteb +QCAPE*q - gamma*(tht1+gamma2*tht2)),0);
CAPEc = max(CAPEbar +R*(thteb - gamma*(tht1 + gamma2p*tht2) ),0);
CAPEd=max(0,(1/tau_conv)*(a1*thteb+a2*q - a0*(tht1+gamma2*tht2)));
tebmtem =   thteb - q - tspi*(tht1+alpha2*tht2);


Hdp = max( fdeq*sqrt(CAPEbar)/Hm + (1/tau_conv)*(a1*thteb+a2*q - a0*(tht1+gamma2*tht2)),0);
Hd  = fd.*Hdp/fdeq;

%Hs =  fs.*sqrt(CAPE)/Hm/4;
if(delay==false)
Hsp = max( fseq*sqrt(CAPEbar)/Hm + (1/tau_conv)*(a1*thteb+a2*q - a0*(tht1+gamma2*tht2)),0);
Hs  = (fs/fseq).*(Hsp)/4;
else
  Hs=sol(:,8);
end
  

Hc =  fc.*sqrt(CAPEc)/Hm*alphac;

P=Hd+xi_c*Hc+xi_s*Hs;
Pl=(Hd+xi_s*Hs+xi_c*Hc);
Plm=mean(Pl);
f_s=mean(xi_s*Hs)/Plm;
f_c=mean(xi_c*Hc)/Plm;
f_d=1-f_c-f_s;



RCEP=(fdeq+xi_c*fceq*alpha_c +xi_s*fseq*alpha_s);
RCEf_s=xi_s*fseq*alpha_s/RCEP;
RCEf_c=xi_c*fceq*alpha_c/RCEP;
RCEf_d=fdeq/RCEP;




figure(1)
subplot1(5,1,'Gap',[0.05 0.035],'FontS',14)


subplot1(1)

h=plot(NDAYS/2+ (0 :Niter)*DT0/24,sol(:,1:2)*alpha_bar);
legend(' $\theta_1$ ',' $\theta_2$ ','Location','NorthEast')
set(h,'LineWidth',2,{'LineStyle'},{'-';'--'})
set(h,{'Color'},{'g';'k'})
ylabel('K','fontsize',14)

ylim([min(min(sol(:,1:2)*alpha_bar)), max(max(sol(:,1:2)*alpha_bar))])
hold on





subplot1(2)
h=plot(NDAYS/2+ (0 :Niter)*DT0/24,tht1*alpha_bar, NDAYS/2+ (0 :Niter)*DT0/24,thteb*alpha_bar,NDAYS/2+ (0 :Niter)*DT0/24,q*alpha_bar);
legend(' $\theta_{1}$',' $\theta_{eb}$','q','Location','NorthEast')
set(h,'LineWidth',2,{'LineStyle'},{'-';'-';'--'})
set(h,{'Color'},{'g';'k';'k'})
ylabel('K','fontsize',14)

ylim([min(min(min(sol(:,3:4)*alpha_bar)),min(tebmtem*alpha_bar)), max(max(max(sol(:,3:4)*alpha_bar)),max(tebmtem*alpha_bar))])
hold on






subplot1(3)

h2=plot(NDAYS/2+ (0 :Niter)*DT0/24,sol(:,5:7));
set(h2,'LineWidth',2)
legend(' $\sigma_c$','$\sigma_d$ ',' $\sigma_s$ ','Location','NorthEast')
ylabel('','fontsize',14)
ylim([min(min(sol(:,5:7))), max(max(sol(:,5:7)))]);
tend =NDAYS/2+ Niter*DT0/24;
set(h2,'LineWidth',2,{'LineStyle'},{'-';'-';'--'})
set(h2,{'Color'},{'g';'k';'k'})

hold on

 



% subplot1(3)
% 
% h2=plot(NDAYS/2+ (0 :Niter)*DT0/24,CAPE*50^2,NDAYS/2+ (0 :Niter)*DT0/24,CAPEc*50^2);
% set(h2,'LineWidth',2)
% legend('CAPE ','CAPEl ',0)
% ylabel('J/kg','fontsize',14)
% hold on
% tend = Niter*DT0/24;
% plot([0 tend],[CAPEbar CAPEbar]*50^2,'b--')
% plot([0 tend],[fseq, fseq],'r--')
% hold off
% 




subplot1(4)

h2=plot(NDAYS/2+ (0 :Niter)*DT0/24,Hc*15.33/(8.33/24),NDAYS/2+ (0 :Niter)*DT0/24,Hd*15.33/(8.33/24),NDAYS/2+ (0 :Niter)*DT0/24,Hs*15.33/(8.33/24));
set(h2,'LineWidth',2)

ylim([0, max(max(max(Hc,Hs),max(Hd)))]*15.33/(8.33/24));



legend(' $H_c$ ',' $H_d$ ',' $H_s$ ','Location','NorthEast')
ylabel('K/Day','fontsize',14)
hold on
tend = Niter*DT0/24;
% plot([0 tend],[fceq fceq],'b--')
% plot([0 tend],[fdeq, fdeq],'g--')
% plot([0 tend],[fseq, fseq],'r--')
%xlabel('days','fontsize',14)
hold on


set(h2,'LineWidth',2,{'LineStyle'},{'-';'-';'--'})
set(h2,{'Color'},{'g';'k';'k'})

title([strcat('RCE:  $f_c$ =',num2str(round(RCEf_c*100)/100))...
    ,',   f_s =',num2str(round(RCEf_s*100)/100),',   $f_d/f_s$ =',num2str(round(RCEf_d/RCEf_s*100)/100),... 
    ', ACTUAL:  $f_c$ =',num2str(round(f_c*100)/100)...
    ,',   $f_s$ =',num2str(round(f_s*100)/100),',   $f_d/f_s$ =',num2str(round(f_d*100/f_s)/100)]);









subplot1(5)
tebmtem = tebmtembar + thteb - q - tspi*(tht1+alpha2*tht2);

h2=plot(NDAYS/2+ (0 :Niter)*DT0/24,max(0,1-exp(-CAPEc/CAPE0)),'k--',NDAYS/2+ (0 :Niter)*DT0/24,max(0,1-exp(-CAPE/CAPE0))...
    ,'k-',NDAYS/2+ (0 :Niter)*DT0/24,max(0,exp(-max(tebmtem,0)/MOIST0)) ,'ro');
set(h2,'LineWidth',2)



ylim=([0, 1]);
legend(' $\Gamma(C_l)$ ','$\Gamma(C)$','$1-\Gamma(D)$','Location','NorthEast')
ylabel('nd','fontsize',14)
xlabel('days','fontsize',14)

set(h2,'LineWidth',2,{'LineStyle'},{'-';'-';'--'})
set(h2,{'Color'},{'g';'k';'k'})
hold on



[peakLoc]=peakfinder(Pl,max(P)/3);

Period=peakLoc(2:end)-peakLoc(1:(end-1));

title([strcat('Mean Precip=',num2str(round((15.33/(8.33/24))*mean(Pl)*100)/100))...
    ,',  Period =',num2str(round((DT0)*mean(Period)*100)/100),'Hours,  	St.Dv. Precip=',num2str(round((15.33/(8.33/24))*std(P)*100)/100),'K/D,  St.Dv. Period=',num2str(round((DT0)*std(Period)*100)/100),'H']);




oldSettings = fillPage(gcf, 'margins', -[.5 2 1.5 2]/3);
print(gcf, '-dpdf', '-r300', 'McLong.pdf')
closeup=input('enterclouseup');
for J=1:6
    subplot1(J)
    xlim(closeup)
    alim
end

oldSettings = fillPage(gcf, 'margins', -[.5 2 1.5 2]/3);
print(gcf, '-dpdf', '-r300', 'Mcshort.pdf')

if(stoch)
    WLINE=2;
else
    WLINE=6;
end
figure(221)
subplot1(2,1,'Gap',[0.05 0.035],'FontS',14)


N=max(size(P));
Timeint=NDAYS/2+ (0 :Niter)*DT0/24;
twk=[Timeint(1),Timeint(end)];
Tend=Timeint(end)-Timeint(1);
%f = sin(2*pi*Timeint)+sin(2*pi*Timeint/5)+randn(size(Timeint)); %%define function, 10 Hz sine wave
f=thteb;
f=f-mean(f);


p = abs(fft(f))/(N/2); %% absolute value of the fft
p = p(1:N/2).^2; %% take the power of positve freq. half
freq = [0:N/2-1]/Tend; %% find the corresponding frequency in Hz


subplot1(1)
plot(freq,p); %% plot on s
%axis tight
axis([0 2.25 0 max(p)])
set(gca,'XTick',[0:.25:2.5])
grid on


f=P;
f=f-mean(f);


p = abs(fft(f))/(N/2); %% absolute value of the fft
p = p(1:N/2).^2; %% take the power of positve freq. half
freq = [0:N/2-1]/Tend; %% find the corresponding frequency in Hz

legend('$\theta_{eb}$. ','Location','NorthEast')
xlabel('frequency 1/days','fontsize',14)

subplot1(2)
plot(freq,p); %% plot on s
%axis tight
axis([0 2.5 0 max(p)])
set(gca,'XTick',[0:.25:2.5])
grid on

legend('Precip ','Location','NorthEast')
xlabel('frequency 1/days','fontsize',14)




%axis([0 2.5 0 1]); %% zoom in
%yaxis tight



