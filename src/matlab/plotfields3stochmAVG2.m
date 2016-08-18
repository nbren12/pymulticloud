
varv1=max(50*var(pv1'));
varv2=max(50*var(pv2'));
vart1=max(15*var(pt1'));
vart2=max(15*var(pt2'));
varteb=max(15*var(pteb'));
varq=max(15*var(pq'));
varHc=max(15/(8.33/24)*var(phc'));
varHd=max(15/(8.33/24)*var(phd'));
varHs=max(15/(8.33/24)*var(phs'));
varFc=max(var(100*pfc'));
varFd=max(var(100*pfd'));


figure (7)
%imagesc
subplot1(4,2,'Gap',[0.0001 0.04],'FontS',14)
subplot1(1)
%contourf(x,time,pl')
contourf(x,time,50*(pv1'-j*mean(pv1')),5)
set(gca,'YDir','normal')
colorbar
ylabel('days')
%title(' v_1  (m/s)')
title([strcat('v_1 m/s | max var. =',num2str(round(varv1*100)/100))])



subplot1(2)
%contourf(x,time,pl')
contourf(x,time,50*(pv2'-j*mean(pv2')),5)
set(gca,'YDir','normal')
colorbar
title([strcat('v_2 m/s | max var. =',num2str(round(varv2*100)/100))])

subplot1(3)
%contourf(x,time,pl')
contourf(x,time,15*(pt1'-j*mean(pt1')),5)
set(gca,'YDir','normal')
colorbar
ylabel('days')
title([strcat('\theta_1 K | max var. =',num2str(round(vart1*100)/100))])

subplot1(4)
%contourf(x,time,pl')
contourf(x,time,15*(pt2'-j*mean(pt2')),5)
set(gca,'YDir','normal')
colorbar
title([strcat('\theta_2 K | max var. =',num2str(round(vart2*100)/100))])


subplot1(5)
%contourf(x,time,pl')
contourf(x,time,15*(pteb'-j*mean(pteb')),5)
set(gca,'YDir','normal')
colorbar
ylabel('days')
title([strcat('\theta_{eb} K | max var. =',num2str(round(varteb*100)/100))])

subplot1(6)
%contourf(x,time,pl')
contourf(x,time,15*(pq'-j*mean(pq')),5)
set(gca,'YDir','normal')
colorbar
title([strcat('q K | max var. =',num2str(round(varq*100)/100))])


subplot1(7)
%contourf(x,time,pl')
contourf(x,time,(100*pfd'-j*mean(pfd')),5,'LineStyle','none')
set(gca,'YDir','normal')
colorbar
xlabel('(x1000)KM ')
ylabel('days')
title([strcat('\sigma_d % | max var. =',num2str(round(varFd*10000)/10000))])

subplot1(8)
%contourf(x,time,pl')
contourf(x,time,(100*pfc'-j*mean(pfc')),5,'LineStyle','none')
set(gca,'YDir','normal')
colorbar
xlabel('(x1000)KM ')
title([strcat('\sigma_c % | max var. =',num2str(round(varFc*10000)/10000))])



oldSettings = fillpage(gcf, 'margins', -[.5 2 1.5 2]/3);
print(gcf, '-dpdf', '-r300', 'ALLmAVG.pdf')




