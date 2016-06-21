function  [pi_eq,r01,r02,r10,r20,r23,r12,r30]=  equilibriumdistribution2fort(C,Cl,D);
global tau01
global tau02
global tau10
global tau12
global tau20
global tau30
global tau23
global r23value
global timeset
global Tau 


    




    r01 = gammabb(Cl)*gammabb(D)/tau01;
    r02 = gammabb(C)*(1.d0-gammabb(D))/tau02;
    r10 = gammabb(D)/tau10;
    r12 = gammabb(C)*(1.d0-gammabb(D))/tau12;
    r20 = (1-gammabb(C))/tau20;
    
     r30 = 1.0/tau30;
     
    switch(r23value)
         case(1)
           r23 = 1./tau23;
         case(2)
           r23 = gammabb(sqrt(C))/tau23;
        case(3)
           r23 = gammabb(C)/tau23;
     end
    
    
    x(1) = 1;
     if(r01==0)
         x(2) = 0;
     else
         x(2) = r01/(r10+r12);
     end
     
     
     
    if(r12*r01==0)
        x(3) = 1/(r20+r23)*(r02);
    else
        x(3) = 1/(r20+r23)*(r02+(r12*r01)/(r10+r12));
    end
    
    
    x(4)=x(3)*r23/r30;
    
    pi_eq=x/sum(x); 
    
    