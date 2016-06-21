

if(avgHD)
v1=50*AVGpv1;
v2=50*AVGpv2;
t1=15.3*AVGpt1;
t2=15.3*AVGpt2;
teb=15.3*AVGpteb;
q=15.3*AVGpq;
Hc=15/(8.33/24)*AVGphc;
Hd=15/(8.33/24)*AVGphd;
Fc=AVGpfc;
Fs=AVGpfs;
Fd=AVGpfd;
Hs=15/(8.33/24)*AVGphs;
Precip=Hd+xi_s*Hs+Hc*xi_c;


else
    
v1=50*mean(pv1');
v2=50*mean(pv2');
t1=15*mean(pt1');
t2=15*mean(pt2');
teb=15*mean(pteb');
q=15*mean(pq');
Hc=15/(8.33/24)*mean(phc');
Hd=15/(8.33/24)*mean(phd');
Fc=mean(pfc');
Fs=mean(pfs');
Fd=mean(pfd');
Hs=15/(8.33/24)*mean(phs');    
end

ebmem=(teb-(q+(2*sqrt(2)/pi)*(t1+.1*t2)));
v1pol=spline(x,v1);
v2pol=spline(x,v2);
v1p=fnder(v1pol);
v2p=fnder(v2pol);
v1p=ppval(v1p,x);
v2p=ppval(v2p,x);





xavg=max(size(AVGpv1));
xavg=linspace(0,40,xavg);

figure (4)
subplot1(6,1,'Gap',[0.05 0.01],'FontS',14)
subplot1(1)
%contourf(x,time,pl')
plot(x,v1,x,v2)
set(gca,'YDir','normal')
ylabel('m/s')
hl=legend('$v_1$','$v_2$');
     set(hl, 'Color', 'none','FontSize',14)
ylim([min(min(v1,v2)), max(max(v1,v2))])
 grid on


subplot1(2)
%contourf(x,time,pl')
plot(x,t1,x,t2)
set(gca,'YDir','normal')
ylabel('K')
hl=legend('$\theta_1$','$\theta_2$');
     set(hl, 'Color', 'none','FontSize',14)
     grid on
ylim([min(min(t1,t2)), max(max(t1,t2))])


subplot1(3)
%contourf(x,time,pl')
plot(x,teb,x,q,x,ebmem)
set(gca,'YDir','normal')
ylabel('K')


hl=legend('$\theta_{eb}$','q','$\theta_{eb}-\theta_{em}$');
set(hl, 'Color', 'none','FontSize',14)

%  set(hl, 'Color', 'none','FontSize',12)
 grid on
ylim([min(min(teb,min(q,ebmem))), max(max(teb,max(q,ebmem)))])


subplot1(4)
%contourf(x,time,pl')
plot(x,Hd,x,Hc,x,Hs,x,Precip)
temp1=15/(8.33/24)*max(max(max(mean(phd'),max(mean(phc')))),mean(phs'));
set(gca,'YDir','normal')
ylabel('K/day')
%hl=legend([strcat('$H_d$| mean=',num2str(round(mean(Hd)*10)/10))],[strcat('$H_c$| mean=',num2str(round(mean(Hc)*10)/10))],[strcat('$H_s$| mean=',num2str(round(mean(Hs)*10)/10))]);
hl=legend('$H_d$','$H_c$','$H_s$','Precip' );
set(hl, 'Color', 'none','FontSize',14)
ylim([0, max(max(max(Hd,Hc)),max(Hs))])

grid on



subplot1(5)
%contourf(x,time,pl')
plot(x,Fc,x,Fd,x,Fs)
temp1=15/(8.33/24)*max(max(mean(pfd'),max(mean(pfc'))));
set(gca,'YDir','normal')
ylabel('')
hl=legend('$\sigma_c$','$\sigma_s$','$\sigma_d$')
set(hl, 'Color', 'none','FontSize',14)
hold on
onz=ones(size(Fd)) ;
  plot(x,onz*fdeq,x,onz*fceq,x,onz*fseq)
ylim([0, max(max(max(Fd,Fc)),max(Fs))])

grid on









z=0:.25:16;
H=zeros(max(size(z)),max(size(x)));
Theta=H;
U=H;
W=H;




for n=1:(max(size(z)))
    
H(n,:)=(Hd-mean(Hd))*F1(z(n))   +   (Hc-mean(Hc))*(F2(z(n))-delc*F2(16-z(n)))+(Hs-mean(Hs))*(F2(16-z(n)) -dels*F2(z(n)) );
Theta(n,:)=15*z(n)*pi/16+sqrt(2)*t1*sin(z(n)*pi/16)+sqrt(2)*2*t2*sin(2*z(n)*pi/16);
%Theta(n,:)=t1*F1(z(n))+2*t2*F2(2*z(n));
%Theta(n,:)=t1*sin(z(n)*pi/16)+2*t2*sin(2*z(n)*pi/16);
Theta(n,:)=sqrt(2)*(t1-mean(t1))*sin(z(n)*pi/16)+sqrt(2)*2*(t2-mean(t2))*sin(2*z(n)*pi/16);

U(n,:)=v1*cos(z(n)*pi/16)+v2*cos(2*z(n)*pi/16);
W(n,:)=-(v1p*sin(z(n)*pi/16)+(1/2)*v2p*sin(2*z(n)*pi/16))*16/pi;
end


zc=z(1:5:end);
xc=x(1:40*xgrid/1040:end);
Uc1=U(1:5:end,:);
Uc=Uc1(:,1:40*xgrid/1040:end);
Wc1=W(1:5:end,:);
Wc=Wc1(:,1:40*xgrid/1040:end);

Wmax=int2str(max(max(abs(W)))/10); %cm/s
Umax=int2str(max(max(abs(U))));  %m/s


subplot1(6)

tmax=max(max(Theta));
tmin= min(min(Theta));
v=tmin:1:tmax;
contour(x,z,Theta,5)
hold on
ylabel('height (km)')
xlabel('(x1000)KM ')
quiver(xc,zc,Uc,Wc)
colorbar('east')

hl=legend([strcat('$\Theta | max(U,W)$ =(',Umax,'m/s',',',Wmax,'cm/s',')    ')],'Location','NorthWest');
     set(hl, 'Color', 'none','FontSize',14)

% 
% subplot1(7)
% tmax=max(max(H));
% tmin= min(min(H));
% v=tmin:2:tmax;
% contour(x,z,H,5)
% hold on
% quiver(xc,zc,Uc,Wc)
% legend('H ','Location','NorthWest')
% ylabel('height (km)')
% xlabel('X (1000km) ')
% colorbar('east')


oldSettings = fillPage(gcf, 'margins', -[.5 2 1.5 2]/3);
print(gcf, '-dpdf', '-r300', 'AVGsnapHD.pdf')
 saveas(gcf, 'AVGsnapHD', 'fig')

