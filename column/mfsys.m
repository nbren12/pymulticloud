function dy = mfsys(t,y)%,r01,r02,r10,r20,r23,r12,r30)
% function to be integrated
global Gr01 Gr02 Gr10 Gr20 Gr23 Gr12 Gr30

nc=1-y(1)-y(2)-y(3);
dy = zeros(3,1);
dy(1) = nc*Gr01-y(1)*(Gr10+Gr12);
dy(2) = nc*Gr02+y(1)*(Gr12)-y(2)*(Gr20+Gr23);
dy(3) = y(2)*Gr23-y(3)*Gr30;
