 function  setRc(C,D)
 
global Rc Tau

global Rc Tau xi_s xi_c alpha_s alpha_c
global r23value
global timeset






Rc(1,2)=G7(C)*G7(D)/Tau(1,2);
Rc(2,1)=G7(D)/Tau(2,1);
Rc(2,3)=G7(C)*(1-G7(D))/Tau(2,3);
Rc(1,3)=G7(C)*(1-G7(D))/Tau(1,3);
Rc(3,4)=1/Tau(3,4);
Rc(4,1)=1/Tau(4,1);
Rc(3,1)=(1-G7(C))/Tau(3,1);

if (r23value==2)
    Rc(3,4)=G7(sqrt(C))/Tau(3,4);

end






