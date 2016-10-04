
varv1=max(50*var(pv1'));
varv2=max(50*var(pv2'));
vart1=max(15*var(pt1'));
vart2=max(15*var(pt2'));
varteb=max(15*var(pteb'));
varq=max(15*var(pq'));
varHc=max(15/(8.33/24)*var(phc'));
varHd=max(15/(8.33/24)*var(phd'));
varHs=max(15/(8.33/24)*var(phs'));
varFc=max(var(pfc'));
varFd=max(var(pfd'));
mono=true;
maxcutoff=30;

J=size(pv1);
J=J(2);
J=ones(1,J);

V1=(v1*J)';
V2=(v2*J)';

T1=(t1*J)';
T2=(t2*J)';

Teb=(teb*J)';
Q=(q*J)';

HC=(Hc*J)';
HD=(Hd*J)';
HS=(Hs*J)';


FC=(Fc*J)';
FD=(Fd*J)';
FS=(Fs*J)';



figure (23)
%imagesc
subplot1(6,2,'Gap',[0.01 0.04],'FontS',14)
subplot1(1)
%contourf(x,time,pl')
contourf(x,time,50*(pv1')-V1,5,'LineStyle','none')
set(gca,'YDir','normal')
colorbar

if(mono)
colormap(bone)
colormap(flipud(colormap))
end
ylabel('days')
%title(' v_1  (m/s)')
title([strcat(' $v_1$  m/s | max var. =',num2str(round(varv1*100)/100))])



subplot1(2)
%contourf(x,time,pl')
contourf(x,time,50*(pv2')-V2,5,'LineStyle','none')
set(gca,'YDir','normal')
colorbar

if(mono)
colormap(bone)
colormap(flipud(colormap))
end
title([strcat(' $v_2$  m/s | max var. =',num2str(round(varv2*100)/100))])

subplot1(3)
%contourf(x,time,pl')
contourf(x,time,15*(pt1')-T1,5,'LineStyle','none')
set(gca,'YDir','normal')
colorbar

if(mono)
colormap(bone)
colormap(flipud(colormap))
end
ylabel('days')
title([strcat(' $\theta_1$  K | max var. =',num2str(round(vart1*100)/100))])

subplot1(4)
%contourf(x,time,pl')
contourf(x,time,15*(pt2')-T2,5,'LineStyle','none')
set(gca,'YDir','normal')
colorbar

if(mono)
colormap(bone)
colormap(flipud(colormap))
end
title([strcat(' $\theta_2$  K | max var. =',num2str(round(vart2*100)/100))])


subplot1(5)
%contourf(x,time,pl')
contourf(x,time,15*(pteb')-Teb,5,'LineStyle','none')
set(gca,'YDir','normal')
colorbar

if(mono)
colormap(bone)
colormap(flipud(colormap))
end
ylabel('days')
title([strcat(' $\theta_{eb}$  K | max var. =',num2str(round(varteb*100)/100))])

subplot1(6)
%contourf(x,time,pl')
contourf(x,time,15*(pq')-Q,5,'LineStyle','none')
set(gca,'YDir','normal')
colorbar

if(mono)
colormap(bone)
colormap(flipud(colormap))
end
title([strcat('q K | max var. =',num2str(round(varq*100)/100))])

%contourf(x,time,pl')


subplot1(8)
%contourf(x,time,pl')
contourf(x,time,(pfd')-FD,20,'LineStyle','none')
set(gca,'YDir','normal')
colorbar
ylabel('days')
title([strcat(' $\sigma_d$   | max var. =',num2str(round(varFd*100)/100))])

subplot1(7)
%contourf(x,time,pl')
contourf(x,time,(pfc')-FC,20,'LineStyle','none')
set(gca,'YDir','normal')
colorbar
title([strcat(' $\sigma_c$   | max var. =',num2str(round(varFc*100)/100))])


subplot1(9)
%contourf(x,time,pl')
contourf(x,time,min(15/(8.33/24)*(phc')-HC,maxcutoff),20,'LineStyle','none')
set(gca,'YDir','normal')
colorbar
title([strcat('$H_c$ K/day | max var. =',num2str(round(varHc*100)/100))])
%
if(mono)
colormap(bone)
colormap(flipud(colormap))
end
xlabel('KM ')

subplot1(10)
contourf(x,time,min(15/(8.33/24)*(phd')-HD,maxcutoff),20,'LineStyle','none')
set(gca,'YDir','normal')
colorbar
ylabel('days')
title([strcat('$H_d$ K/day | max var. =',num2str(round(varHd*100)/100))])
%
if(mono)
colormap(bone)
colormap(flipud(colormap))
end
xlabel('KM ')




subplot1(11)
%contourf(x,time,pl')
contourf(x,time,pfs'-FS,20,'LineStyle','none')
set(gca,'YDir','normal')
colorbar
title([strcat('$\sigma_s$ K/day | max var. =',num2str(round(varHc*100)/100))])
%
if(mono)
colormap(bone)
colormap(flipud(colormap))
end
xlabel('KM ')

subplot1(12)
contourf(x,time,min(15/(8.33/24)*(phs')-HS,maxcutoff),20,'LineStyle','none')
set(gca,'YDir','normal')
colorbar
ylabel('days')
title([strcat('$H_s$ K/day | max var. =',num2str(round(varHd*100)/100))])
%
if(mono)
colormap(bone)
colormap(flipud(colormap))
end
xlabel('KM ')








oldSettings = fillPage(gcf, 'margins', -[.5 1 1.5 2]/3);
print(gcf, '-dpdf', '-r300', 'AllsnapMAVG.pdf')




