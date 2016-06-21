function     [fc,fd,fs,count] = birthdeathexact2fort(fc,fd,fs,C,Cl,D,DT, n)


% Number of sites = n^2
nsite = n^2;

% Conditional rates
[pi_eq,r01,r02,r10,r20,r23,r12,r30]=  equilibriumdistribution2fort(C,Cl,D);
% Define cloud state
cloud = round(nsite*[fc,fd,fs]);
%min(cloud)

% define the seven possible transitions
diffv(1,1:3) = [1,0,0];
diffv(2,1:3) = [-1,0,0];
diffv(3,1:3) = [-1,1,0];
diffv(4,1:3) = [0,1,0];
diffv(5,1:3) = [0,-1,0];
diffv(6,1:3) = [0,-1,1];
diffv(7,1:3) = [0,0,-1];

% We will use Gilespie's algorithm and iterate until we have exceeded the
% time step, DT

time = 0;
count = 0;
while time < DT
    count = count+1;
    % Calculate the number of clearsky elements
    clearsky = nsite - sum(cloud);


    % Absolute rates 
    R01 = r01*clearsky;
    R10 = r10*cloud(1);
    R12 = r12*cloud(1);
    R02 = r02*clearsky;
    R20 = r20*cloud(2);
    R23 = r23*cloud(2);
    R30 = r30*cloud(3);
    
    Rmat = [R01,  R10,  R12,  R02,  R20,  R23,  R30];
    
    Rsum = cumsum(Rmat);
    
    time;
    if(Rsum(7) < 0)
        'RSUM ERROR'
            Rmat
             Rsum
             count
             time
             cloud
    
        pause
    elseif(Rsum(7) == 0)
        t = 2*DT;
    elseif(Rsum(7) > 0)
        t = -1*log(rand(1,1))/Rsum(7);

    end    

    
    % Calculate the time until the next transition

    

    time = time+t;
    
%%reject transition if time exceeds DT 
    if(time>DT)
        break
    end
    
    % Which transition occurs?
    test = rand(1,1);
    Rsum = Rsum/Rsum(7); 
    Rsum = Rsum - test;


    
    tmp =  find(sign(Rsum)==1);
    
    ind = min(tmp);
    diff = sum(diffv(ind,:),1);


    
   


    % New cloud array
    cloud = cloud + diff;
    

    
    
end


fc = cloud(1)/nsite;
fd = cloud(2)/nsite;
fs = cloud(3)/nsite;

