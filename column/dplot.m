
WLINE=2;
sold=zeros(288,8);

CAPEd=zeros(288,1);
CAPEcd=zeros(288,1);
Hcd=zeros(288,1);
Hdd=zeros(288,1);
dday=linspace(0,24,288);
for I=1:288
    sold(I,:)=mean(sol(I:288:end,:));
CAPEd(I)=mean(CAPE(I:288:end,:));
CAPEcd(I)=mean(CAPEc(I:288:end,:));
Hcd(I)=mean(Hc(I:288:end,:));
Hdd(I)=mean(Hd(I:288:end,:));
end
tht1=sold(:,1);
tht2=sold(:,2);
thteb=sold(:,3);
q=sold(:,4);
fc=sold(:,5);
fd=sold(:,6);
if(delay)
fs=sold(:,7);
else
end

tebmtem =   thteb - q - tspi*(tht1+alpha2*tht2);




%Hs =  fs.*sqrt(CAPE)/Hm/4;
if(delay==false)
Hsp = max( fseq*sqrt(CAPEbar)/Hm + (1/tau_conv)*(a1*thteb+a2*q - a0*(tht1+gamma2*tht2)),0);
Hs  = (fs/fseq).*(Hsp)/4;
else
  Hs=sold(:,8);
end
  


P=Hdd+xi_c*Hcd+xi_s*Hs;
Pl=(Hdd+xi_s*Hs+xi_c*Hcd);
Plm=mean(Pl);
f_s=mean(xi_s*Hs)/Plm;
f_c=mean(xi_c*Hc)/Plm;
f_d=1-f_c-f_s;



RCEP=(fdeq+xi_c*fceq*alpha_c +xi_s*fseq*alpha_s);
RCEf_s=xi_s*fseq*alpha_s/RCEP;
RCEf_c=xi_c*fceq*alpha_c/RCEP;
RCEf_d=fdeq/RCEP;



figure(11)
hold on
subplot1(5,1,'Gap',[0.05 0.035],'FontS',14)


subplot1(1)

h=plot(dday,sold(:,1:2)*alpha_bar);
legend(' $\theta_1$ ',' $\theta_2$ ','Location','NorthEast')
set(h,'LineWidth',WLINE,{'LineStyle'},{'-';'--'})
set(h,{'Color'},{'g';'k'})
ylabel('K','fontsize',14)


hold on
axis tight



subplot1(2)
h=plot(dday,thteb*alpha_bar,dday,q*alpha_bar);
legend(' $\theta_{eb}$ ','q','Location','NorthEast')
set(h,'LineWidth',WLINE,{'LineStyle'},{'-';'--'})
set(h,{'Color'},{'g';'k'})
ylabel('K','fontsize',14)

hold on
axis tight


xlim([0,24])




subplot1(3)

h2=plot(dday,sold(:,5:7));
set(h2,'LineWidth',2)
legend(' $\sigma_c$ ',' $\sigma_d$ ',' $\sigma_s$ ','Location','NorthEast')
ylabel('','fontsize',14)
axis tight
tend =NDAYS/2+ Niter*DT0/24;
set(h2,'LineWidth',WLINE,{'LineStyle'},{'-';'-';'--'})
set(h2,{'Color'},{'g';'k';'k'})

hold on

 
xlim([0,24])



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

h2=plot(dday,Hcd*15.33/(8.33/24),dday,Hdd*15.33/(8.33/24),dday,Hs*15.33/(8.33/24));
set(h2,'LineWidth',2)




legend(' $H_c$ ',' $H_d$ ',' $H_s$ ','Location','NorthEast')
ylabel('K/Day','fontsize',14)
hold on
tend = Niter*DT0/24;
% plot([0 tend],[fceq fceq],'b--')
% plot([0 tend],[fdeq, fdeq],'g--')
% plot([0 tend],[fseq, fseq],'r--')
%xlabel('days','fontsize',14)
hold on


set(h2,'LineWidth',WLINE,{'LineStyle'},{'-';'-';'--'})
set(h2,{'Color'},{'g';'k';'k'})

title('diurnal deviations from the mean');

axis tight

xlim([0,24])






subplot1(5)
tebmtem = tebmtembar + thteb - q - tspi*(tht1+alpha2*tht2);


h2=plot(dday,max(0,1-exp(-CAPEcd/CAPE0)),'k--',dday,max(0,1-exp(-CAPEd/CAPE0))...
    ,'k-',dday,max(0,exp(-max(tebmtem,0)/MOIST0)) ,'ro');


legend(' $\Gamma(C_l)$ ','$\Gamma(C)$','$1-\Gamma(D)$','Location','NorthEast')
ylabel('nd','fontsize',14)
xlabel('Hours','fontsize',14)

set(h2,'LineWidth',WLINE,{'LineStyle'},{'-';'-';'--'})
set(h2,{'Color'},{'g';'k';'k'})
hold on
axis tight

xlim([0,24])


[peakLoc]=peakfinder(Pl,max(P)/3);

Period=peakLoc(2:end)-peakLoc(1:(end-1));





oldSettings = fillPage(gcf, 'margins', -[.5 2 1.5 2]/3);
print(gcf, '-dpdf', '-r300', 'dplot.pdf')
