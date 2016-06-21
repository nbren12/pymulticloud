function     [fc,fd,fs,count] = mf(fc,fd,fs,C,Cl,D,DT,n,T)



% Conditional rates
global Gr01 Gr02 Gr10 Gr20 Gr23 Gr12 Gr30
[pi_eq,r01,r02,r10,r20,r23,r12,r30]=  equilibriumdistribution2fort(C,Cl,D);

Gr01=r01;
Gr02=r02;
Gr10=r10;
Gr20=r20;
Gr23=r23;
Gr12=r12;
Gr30=r30;



[T,yl]=ode23(@mfsys,[0,DT],[fc,fd,fs]);
fc=yl(end,1);
fd=yl(end,2);
fs=yl(end,3);