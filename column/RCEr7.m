 function [F] = RCEr7(x,QR01,tt,Hm,CAPE0)
 
global Rc Tau xi_s xi_c alpha_s alpha_c
global r23value
global timeset

sigC=x(1);
sigD=x(2);
sigS=x(3);
CAPEb=x(4);

C=CAPEb/CAPE0;
D = tt;


Rc(1,2)=G7(C)*G7(D)/Tau(1,2);
Rc(2,1)=G7(D)/Tau(2,1);
Rc(2,3)=G7(C)*(1-G7(D))/Tau(2,3);
Rc(1,3)=G7(C)*(1-G7(D))/Tau(1,3);
Rc(3,4)=1/Tau(3,4);
Rc(4,1)=1/Tau(4,1);
Rc(3,1)=(1-G7(C))/Tau(3,1);
Rc

if (r23value==2)
    Rc(3,4)=G7(sqrt(C))/Tau(3,4);

end

F(1)=(1-sigC-sigD-sigS)*Rc(1,2)-sigC*(Rc(2,1)+Rc(2,3));
F(2)=(1-sigC-sigD-sigS)*Rc(1,3)-sigD*(Rc(3,1)+Rc(3,4))+sigC*(Rc(2,3));
F(3)=sigD*Rc(3,4)-sigS*Rc(4,1);
Q=sqrt(CAPEb)/Hm;
F(4)=QR01-sigD*Q-xi_s*alpha_s*sigS*Q-xi_c*alpha_c*sigC*Q;

F


