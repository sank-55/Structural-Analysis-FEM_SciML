%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   ONLY NEED TO TUNE THE EKF for the states  PROPERLY EVERYTHING IS FINE 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  This one with most effective observable points %%%%%%%%%%%%

clear all;
close all; clc;
p=1000;
ro=2700; % ( in kg/m3)
yo =70*10^9;%youngs modulus
neu=0.3;
h=0.001;%thickness
a=0.6; %( m)
b=0.4; % (m)
divx = 3;  %input(' enter the no of divisions on the breadth ');
% ex=length of element in x direction
ex=a/(divx);
divy=3;  %input(' enter the no of divisions on the length ');
% ey=length of element in y direction
ey=b/(divy);
%% Natural frequencies for the above 
%% NB: when frequencies are around 8,9 w1 we find scheme 2 is better than the other model
w_1 = 45*pi;   %1*137.2892;
w_6 =  897.6603; 

N = 50000;
deltat = 0.00005;
nObserv= 5 ;       % No of observable zones ( sensor is acting) 

syms x y



Nm=6;
Nn=Nm;
%phi=[sin(pi*x/a)*sin(pi*y/b) sin(3*pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(3*pi*y/b) sin(3*pi*x/a)*sin(3*pi*y/b) sin(4*pi*x/a)*sin(4*pi*y/b)];
phi=[sin(pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(2*pi*y/b) sin(2*pi*x/a)*sin(pi*y/b) sin(2*pi*x/a)*sin(2*pi*y/b) sin(1*pi*x/a)*sin(3*pi*y/b) sin(3*pi*x/a)*sin(1*pi*y/b)];
%phi=[sin(pi*x/a)*sin(pi*y/b)  sin(3*pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(3*pi*y/b)  sin(3*pi*x/a)*sin(3*pi*y/b)  sin(pi*x/a)*sin(5*pi*y/b)  sin(5*pi*x/a)*sin(1*pi*y/b)];

dxxphi=diff(phi,x,2);
dxyphi=diff(diff(phi,x),y);
dyyphi=diff(phi,y,2);
D=(yo*h^3)/(12*(1-neu^2));

Ndam=6;  % Total no. of damage zones ( each having own fd) not the sensors 
% damdox=[0 a/4;a/4 a/2; a/2 3*a/4;3*a/4 a;0 a/4;a/4 a/2];
% damdoy=[0 b/4;0 b/4;0 b/4;0 b/4;b/4 b/2;b/4 b/2;b/4 b/2];

damdox=[0 a/3; a/3 2*a/3;2*a/3 a;0 a/3; a/3 2*a/3;2*a/3 a];
damdoy=[0 b/2;0 b/2;0 b/2;b/2 b;b/2 b;b/2 b];

% Most obsv Sensor location for 6 sensor placement 
% sensox = [0.2949  0.1627  0.4271  0.1525  0.4475  0.3051];
% sensoy = [ 0.1949 0.1026  0.1026  0.2769  0.2769  0.3179];




% Middle term is added as the 
 % sensox=[a/6 a/2 5*a/6  a/6  5*a/6   a/2 ];
 % sensoy=[b/4 b/4 b/4  3*b/4  3*b/4  b/2 ];


 % Effective for 5 sensor location 
    sensox = [0.2949  0.1627  0.4271  0.1525  0.4475 ];
    sensoy = [ 0.1949 0.1026  0.1026  0.2769  0.2769 ];
%  sensox=[  a/6   5*a/6   1*a/6   5*a/6    a/2   ];
%  sensoy=[   b/4   b/4    3*b/4    3*b/4    b/2 ];


for i=1:Nm
    for j=1:Nm
        Mgm(i,j)=double(ro*h*int(int(phi(i)*phi(j),x,0,a),y,0,b));
        Kgk(i,j)=double(D*int(int(dxxphi(i)*dxxphi(j)+2*dxyphi(i)*dxyphi(j)+dyyphi(i)*dyyphi(j),x,0,a),y,0,b));
    end
    Fgf(i,1)=double(int(int(p*phi(i),x,0,a),y,0,b));
end
   Cgc=0.0003*Mgm+0.0003*Kgk;     

for i=1:Ndam
    for j=1:Nm
        for k=1:Nm
            Kd(i,j,k)=double(D*int(int(dxxphi(j)*dxxphi(k)+2*dxyphi(j)*dxyphi(k)+dyyphi(j)*dyyphi(k),x,damdox(i,1),damdox(i,2)),y,damdoy(i,1),damdoy(i,2)));
        end
    end
end
    

clear x
gamma=1/2;
beta=1/2;

x(:,1)=zeros(Nm,1);
v(:,1)=zeros(Nm,1);
acc(:,1)=zeros(Nm,1);
ac(:,1)=zeros(Nm,1);
Fg(:,1)=zeros(Nm,1);

Kgun=Kgk;  % kGK Is the undamaged one 
Haf=[1 0 0 0];

damf=zeros(Ndam,N);

  for i=200:N
      damf(1,i)=(i-200)^.7*0.3/(N-200)^.7;
  end
% % 
  for i=1200:N
      damf(2,i)=(i-1200)*0.6/(N-1200);
  end
 
  for i=600:N
     damf(3,i)=(i-600)*0.12/(N-600);
 end

% for i=10:N
%     damf(4,i)=((i-10)^0.99)*0.3/(N-10);
% end
% 
for i=30:N
    damf(5,i)=(i-30)*0.45/(N-30);
end


for i=300:N
    damf(6,i)=(i-300)*0.2/(N-300);
end

% % % NON LINEAR DAMAGES 
%   for i=2000:N
%       damf(1,i)=(i-2000)^.7*0.3/(N-2000)^.7;
%   end
% % % 
%   for i=12000:N
%       damf(2,i)=(i-12000)*0.6/(N-12000);
%   end
% %  
%   for i=600:N
%      damf(3,i)=((i-600)^0.8)*0.12/(N-600);
%  end
% 
% for i=100:N
%     damf(4,i)=(i-100)*0.3/(N-100);
% end
% 
% for i=300:N
%     damf(5,i)=(i-300)*0.2/(N-300);
% end
% 
% 
% for i=300:N
%     damf(6,i)=(i-300)*0.2/(N-300);
% end



dim=2*Nm;



for i=1:nObserv
    x1=sensox(i);
    y1=sensoy(i);
   % Phi(i,:)=[sin(pi*x1/a)*sin(pi*y1/b) sin(2*pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(2*pi*y1/b) sin(2*pi*x1/a)*sin(2*pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b) sin(3*pi*x1/a)*sin(3*pi*y1/b)  sin(pi*x1/a)*sin(4*pi*y1/b) sin(2*pi*x1/a)*sin(4*pi*y1/b) sin(3*pi*x1/a)*sin(4*pi*y1/b) sin(4*pi*x1/a)*sin(1*pi*y1/b) sin(4*pi*x1/a)*sin(2*pi*y1/b) sin(4*pi*x1/a)*sin(3*pi*y1/b)];
     %Phi(i,:) = [sin(pi*x1/a)*sin(pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b)  sin(3*pi*x1/a)*sin(3*pi*y1/b)  sin(pi*x1/a)*sin(5*pi*y1/b)  sin(5*pi*x1/a)*sin(1*pi*y1/b)];
     Phi(i,:) = [ sin(pi*x1/a)*sin(pi*y1/b) sin(1*pi*x1/a)*sin(2*pi*y1/b) sin(2*pi*x1/a)*sin(1*pi*y1/b) sin(2*pi*x1/a)*sin(2*pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b)];
end
maxVal = max(Phi, [], 'all');
Phi = Phi./maxVal;

X=zeros(nObserv,Nm);

H=[Phi  X ; X  Phi];
Haf=Phi;


 P=[ones(Nm,1)*10^(-4);ones(Nm,1)*10^(-4)];   % 10^-10  for in 4 sensors
 P=diag(P);
 Pdkf2=P;
%P = 1e-5*((3.8302e-61)^-1)*blkdiag(1.3013e-65,2.6286e-69, 3.4733e-66, 2.8234e-70, 1.6619e-73,6.9362e-69 , 3.2922e-61 , 6.7425e-66, 3.8302e-61, 3.1583e-65, 5.8089e-70, 1.1842e-63);   % for 6 sensor and usual mode shapes
%Pdkf2= 1e-5*((3.8302e-61)^-1)*blkdiag( 1.3047e-65,2.6293e-69, 3.4802e-66, 2.8268e-70, 1.6648e-73,6.9442e-69 , 3.2971e-61 , 6.7527e-66, 3.8364e-61, 3.1618e-65, 5.8148e-70, 1.1854e-63); % for 6 sensor and usual mode shapes

%Pf=eye(Ndam)*10000;  % BEST AT 10000 
Pf = 1e-12*blkdiag(0.3728,0.3122,0.4581,0.3710,0.4272,0.3886);% for 6 sensor scheme 1;
Pf2 = 1e-12*blkdiag(0.3739 ,0.3122,0.7081, 0.6620,  0.6101 , 0.7058); % for 6 sensor scheme 2


% --- 5 sensor based 6 Ndam --------- 
% Pf = 1e-20*blkdiag(0.5549,0.5216, 0.6586 ,0.5991 , 0.5195 ,  0.6397); % for 4 sensor 
% Pf2=1e-20*blkdiag(0.3278, 0.0442, 0.4596 ,0.3721 ,0.4276 ,  0.3899);



%Pf=load('Pf.txt');
 Q=eye(dim)/10^-14;  %  % scheme 1 -> eye(dim)/10^-10 for 6 sensor (unchange modes)
 Qdkf2=eye(dim)/10^-14; % scheme 2 -> eye(dim)/10^-8 for 6 sensor (unchange modes)
%Q = 1e5*P; % for 6 sensor unchange modes
%Qdkf2=1e5*Pdkf2;
%Q(2*Nn+1:dim,2*Nn+1:dim)=eye(dim-(2*Nn));
dimf=1;

%Pekf = 1e-10*eye(2*Nm + Ndam);  %working well in 6s (usual modes) 
%Pekf = 1e6*blkdiag(3.3898e-15,3.8276e-17, 2.4722e-16, 1.3198e-17, 1.8844e-18,2.3570e-17, 3.4714e-08 ,3.8733e-10, 2.5383e-09,  1.3562e-10, 1.8996e-11, 2.4181e-10, 1.9583e-13, 1.3443e-13 ,  2.0059e-13, 1.9658e-13, 1.3482e-13, 2.0093e-13  ); % gives the modes 1e1*multiple
%Qekf = 1e-8*eye(2*Nm + Ndam);
Q_ekf_fd = 1e-8*blkdiag(0.5078,0.5143,0.5107, 0.8089, 0.9138, 0.4110 );      %    0.2663 , 0.2135 , 0.3137, 0.1538 ,0.4348, 0.8606);
P_state_ekf = 1e0*blkdiag( 2.0851e-14,3.4895e-16  , 1.5246e-15 , 4.0725e-17 , 1.1409e-17 ,7.2615e-17, 2.3166e-07 ,3.8771e-09, 2.2563e-09,  4.5250e-10 , 1.2676e-10,  8.0682e-10 );
% Pekf =  blkdiag(1e0*P_state_ekf,Q_ekf_fd);

Pekf = 10^-1*blkdiag(1.3499e-15 , 1.6672e-15 , 9.6627e-17 , 4.6831e-18 ,1.2994e-18 , 1.0910e-18,1.1649e-17 ,1.4995e-08, 1.8522e-10, 1.0734e-09,3.6870e-11 , 1.2121e-11 ,5.2941e-11, 8.4701e-11, 1.9810e-11,2.6663e-11 ,2.6273e-11,8.6801e-10 );

Qekf = Pekf; %1e-4*eye(2*Nm + Ndam); %1e6*blkdiag( 2.0851e-14, 3.4895e-16 , 1.5246e-15, 4.0725e-17, 1.1409e-17 ,7.2615e-17, 2.3166e-07 ,3.8771e-09, 2.5383e-09, 1.6939e-08, 4.5250e-10, 2.4181e-10, 1.9583e-13, 1.3443e-13 ,  2.0059e-13, 1.9658e-13, 1.3482e-13, 2.0093e-13  );

%Qf=    ;
%Qf = 1e-11*blkdiag( 0.2438  ,0.1030 ,0.2463  ,0.2365   , 0.1000, 0.2328);  % working better in 1e-12
Qf = Pf; %1e-12*blkdiag( 0.4793  , 0.3464 ,0.5248, 0.4994 ,0.3564, 0.5304);
Qf2 =Pf2; %    % 1e1*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);

R=1*eye(2*nObserv)*10^(4);  % for state ( disp & vel) of Scheme - 1 for 6 sen ->1*eye(2*nObserv)*10^(-4) works better ( no noise)  
Rf=1*eye(nObserv)*10^(2);
Rf2=1*eye(nObserv)*10^(2);
%Rf = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);

Rdkf2=1*eye(2*nObserv)*10^(4); % for state ( disp & vel) of Scheme - 2
%Rf2 = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);

Rekf = 1e-2*eye(nObserv);

%% nO OF THE sensor reduced Increase the R value

xx=zeros(dim,1);
xx(1,1)=0.0;

Fo(1)=0;

FF(:,1)=zeros(1,Nm);
aa(:,1)=zeros(1,Nm);
vv(:,1)=zeros(1,Nm);

fac=zeros(Ndam,1);
facc=fac;
fd_ex = zeros(Ndam,1);
fd_est_ex = zeros(Ndam,N);
fd_dkf2 = zeros(Ndam,1);
fd_est_dkf2 = zeros(Ndam,N);

xxp=zeros(2*Nm,1);
xxp2=zeros(2*Nm,1);
xxp_exfd=zeros(2*Nm+Ndam,1);  % for the extended kalman filter 
vvp=zeros(Nm,1);
aap=zeros(Nm,1);
accp=zeros(Nm,1);
   acp=zeros(Nm,1);
   acp_2 = zeros(Nm,1);
   xp=zeros(Nm,1);
   vp=zeros(Nm,1);
   xrip=zeros(Nm,1);
   vrip=zeros(Nm,1);
   accpri=zeros(Nm,1);
xx_true_store = zeros(Nm,N);   
xx_store_dkf = zeros(2*Nm,N);
xx_store_dkf2 = zeros(2*Nm,N);
xx_store_ekf = zeros(2*Nm,N);

%% Creating some obs matrix for the different schemes and the states ( sc1 & sc2 & extened) 
Obsv_state_s1 = zeros(2*Nm,1);
Obsv_fd_s2 = zeros(Ndam,1);
Obsv_fd_s1 = zeros(Ndam,1);
Obsv_state_s2 = zeros(2*Nm,1);
Obsv_ekf = zeros(2*Nm + Ndam, 1);

%% Adjusting noise level 
noise_acc = 0.03;
noise_state = 0.02;

for i=2:N
    
     Kgk=Kgun;
%     
    for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgk=Kgk-damf(j,i)*tem;
    end
   
    AA=Mgm+Cgc*deltat*gamma+Kgk*beta*deltat^2;

BB=Cgc*deltat*(1-gamma)+Kgk*(1-2*beta)*deltat^2/2;

CC=Cgc+deltat*Kgk;
AAinv=inv(AA);
    
    t=(i-1)*deltat;
    t;
    Fg(:)=Fgf*sin(w_1*t);
    acc(:)=AAinv*(Fg(:)-[BB]*accp(:)-[CC]*v(:)-Kgk*x(:));
    accp=acc;
    
    for kk=1:Nm
        ac(kk)=acc(kk)+ 0*(rand-.5)/50;
        Accel(i,kk)=acc(kk);
    end
    accri=ac+ 0*(rand-.5)/70;
    accri(1)=ac(1)+ 0*(rand-0.5)/5;
    v(:)=vp+deltat*(1-gamma)*accpri+deltat*gamma*accri;
    x(:)=xp+deltat*vp+(1-2*beta)/2*deltat^2*accpri+beta*deltat^2*accri;
    XX(:,i)=x;
    AC(:,i)=acc;
    ACCRI(:,i)=accri;
    vp=v;
    xp=x;
    
    
    for j=1:Nm
        
        vri(j)=vrip(j)+(acp(j)+ac(j))/2*deltat;
        
        xri(j)=xrip(j)+(vrip(j)+vri(j))/2*deltat;
    end

    xx_true_store(:,i) = xri';
    
    %% Adding noise 
    
     xrri=Phi*(xri + xri*noise_state*rand)';
    vrri=Phi*(vri + vri*noise_state*rand)';
    zf(1:nObserv)=(Phi*(acc + acc*noise_acc*rand))';
    z = [xrri; vrri]';
   



    
    iii=i;
    t=(i-1)*deltat;
    
   
    disp=xxp(1:Nm)+xxp(Nm+1:2*Nm)*deltat;
    vel=xxp(Nm+1:2*Nm)+acp(1:Nm)*deltat;
    fac=facc;
    Kdx=zeros(Nm,Ndam);
  
    for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        xi(:,j)=tem*disp;
        Kdx(:,j)=xi(:,j);
    end
%       F=Fgf*sin(2*pi*t/6);
    F=Fgf*sin(w_1*t);    
    F1=F;
    Mgminv=inv(Mgm);
    
    
    Af=eye(Ndam);
    Hf=Haf*Mgminv*Kdx;
    Cgcx=Mgminv*(Kgun*disp+Cgc*vel);
    HBf=Haf*(Mgminv*F-Cgcx);

    O_s1_fd = obsv(Af,Hf);
    
    % for Scheme 2 

     Kgks2 = Kgun;
     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgks2=Kgks2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
     end

      Kdx2=zeros(Nm,Ndam);
  
    for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        xi(:,j)=tem*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2 );
        Kdx2(:,j)=xi(:,j);
    end
    Hdkf2 = Phi*(inv(Mgm + deltat*0.5*Cgc + deltat^2*0.25*(Kgks2))*(Kdx2)) ;
    a_dkf2 = Hdkf2*fd_dkf2 + Phi*(inv(Mgm + deltat*0.5*Cgc + deltat^2*0.25*(Kgks2))*(F1- Kgun*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2) - Cgc*(xxp2(Nm+1:2*Nm) + deltat^2*0.5*acp_2) ));
    
   O_s2_fd = obsv(Af,Hdkf2);
    % for scheme 1 
    Pnf=Af*Pf*transpose(Af)+Qf;
    Kf=Pnf*transpose(Hf)*inv((Hf*Pnf*transpose(Hf))+Rf);

    %Af=eye(1);
    reba=100;
    fac=facc+Kf*(zf'-Hf*fac(:)-HBf);

     %for Scheme 2
     Pnf2=Af*Pf2*transpose(Af)+Qf2;
    Kdkf2=Pnf2*transpose(Hdkf2)*inv((Hdkf2*Pnf2*transpose(Hdkf2))+Rf2);

    %Af=eye(1);
    reba=100;
    fd_dkf2=fd_dkf2+Kdkf2*(zf'-a_dkf2);






%% 2nd Loop for DEKF 

    Kgk1 = Kgun;
     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgk1=Kgk1-fac(j)*tem;  % here must be fac or facc damf will not be there
     end


    Ac=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgk1  -inv(Mgm)*Cgc];
     Bc=[zeros(Nm,Nm); inv(Mgm)];

 % FOr Scheme 2 

    Kgdkf2 = Kgun;
     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgdkf2=Kgdkf2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
     end


    Adkf2=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgdkf2  -inv(Mgm)*Cgc];
     Bdkf2=[zeros(Nm,Nm); inv(Mgm)];




% -------- Extended kalman filter matrix --------
     S_ekf=zeros(Nm,Ndam);
    for j=1:Ndam
       for k=1:Nm
           for m=1:Nm
               tem(k,m)=Kd(j,k,m);
           end
       end
       S_ekf(:,j) = tem*xxp_exfd(1:Nm);
    end
 
    Aex = [zeros(Nm,Nm) , eye(Nm,Nm) , zeros(Nm,Ndam);
            -inv(Mgm)*Kgk1 , -inv(Mgm)*Cgc , inv(Mgm)*S_ekf ;
            zeros(Ndam,Nm), zeros(Ndam,Nm), eye(Ndam,Ndam)];

    Bex = [zeros(Nm,Nm); inv(Mgm); zeros(Ndam,Nm)];


% Discretization 
% for scheme 1
A=deltat*Ac+eye(2*Nm);
B=deltat*Bc;

% How this is defined
A=zeros(2*Nm,2*Nm);
BB=zeros(2*Nm,2*Nm);
B=zeros(2*Nm,Nm);
for jj=0:10
    A=A+(Ac*deltat)^jj/factorial(jj);
    BB=BB+(Ac*deltat)^jj/factorial(jj+1);
end

B=BB*Bc*deltat;



     %A=exp(Ac*deltat);
     %B=[A-eye(2*Nn)]*inv(Ac)*Bc;
%      A=inv(eye(2*Nm)-Ac);
%      B=inv(eye(2*Nm)-Ac)*Bc;

% for scheme 2 
Ax2=deltat*Adkf2+eye(2*Nm);
Bx2=deltat*Bdkf2;

% How this is defined
Ax2=zeros(2*Nm,2*Nm);
Bx2_=zeros(2*Nm,2*Nm);
Bx2=zeros(2*Nm,Nm);
for jj=0:10
    Ax2=Ax2+(Adkf2*deltat)^jj/factorial(jj);
    Bx2_=Bx2_+(Adkf2*deltat)^jj/factorial(jj+1);
end

Bx2=Bx2_*Bdkf2*deltat;


% ----- For extended kalman filter 
    Aex =  deltat*Aex + eye(2*Nm + Ndam);
    Bex=deltat*Bex;

    Ax=zeros(2*Nm+Ndam,2*Nm+Ndam);
    BBx=zeros(2*Nm+Ndam,2*Nm+Ndam);
    Bx=zeros(2*Nm+Ndam,Nm);
    for jj=0:10
        Ax=Ax+(Aex*deltat)^jj/factorial(jj);
        BBx=BBx+(Aex*deltat)^jj/factorial(jj+1);
    end

    Bx=BBx*Bex*deltat;

    xx_exfd = Ax*xxp_exfd + Bx*F;
    Pekfn = Ax*Pekf*transpose(Ax)+ Qekf; 
     Kgk2 = Kgun;

     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgk2=Kgk2-xx_exfd(2*Nm+j)*tem;  % here must be fac or facc damf will not be there
    end
      S_ekf=zeros(Nm,Ndam);
    for j=1:Ndam
       for k=1:Nm
           for m=1:Nm
               tem(k,m)=Kd(j,k,m);
           end
       end
       S_ekf(:,j) = tem*xx_exfd(1:Nm);
    end
 


    Hekf = [-Phi*inv(Mgm)*Kgk2 , -Phi*inv(Mgm)*Cgc , Phi*inv(Mgm)*S_ekf ];
    z_diff_ekf = zf' - Phi*( Mgm\(F1 - Kgk2*xx_exfd(1:Nm) - Cgc*xx_exfd(Nm+1 :2*Nm) )) ;
    Kekf = Pekfn*transpose(Hekf)*inv((Hekf*Pekfn*transpose(Hekf)) + Rekf);
    xxp_exfd = xx_exfd + Kekf*z_diff_ekf;

    % Storing the Values of the updated ekf 
    xx_store_ekf = xxp_exfd(1:2*Nm); %(  first Nm 'u' then Nm 'vel')
    fd_est_ex(:,i) = xxp_exfd(2*Nm +1: 2*Nm + Ndam);

    O_ex = obsv(Ax,Hekf);

 
   % 0*B*Fg(:,i-1)
   % For Scheme 1 
    xx=A*xxp+B*F;
    Pn=A*P*transpose(A)+Q;    
    K=Pn*transpose(H)*inv((H*Pn*transpose(H))+R);
    xx=xx+K*(z'-H*xx);    
    xx_store_dkf(:,i)=xx;
    vv=(xx(1:Nn)-xxp(1:Nn))/deltat;
    aa=(vv(1:Nn)-vvp(1:Nn))/deltat;
    FF=Mgm*aa(:)+Cgc*vv(:)+Kgk*xx(1:Nn);
    %KK(:,i)=K;
     O_s1_st = obsv(Ac,H);
    
    
    %Fo(i)=Fo(i-1)+deltat/2*(xx(5,i-1)+xx(5,i));
    
    % For Scheme 2 
    xx2=Ax2*xxp2+Bx2*F;
    Pndkf2=Ax2*Pdkf2*transpose(Ax2)+Qdkf2;    
    Kx2=Pndkf2*transpose(H)*inv((H*Pndkf2*transpose(H))+Rdkf2);
    xx2=xx2+Kx2*(z'-H*xx2);   
    xx_store_dkf2(:,i) = xx2;
    fd_est_dkf2(:,i) = fd_dkf2;
    % vv2=(xx2(1:Nn)-xxp2(1:Nn))/deltat;
    % aa2=(vv(1:Nn)-vvp2(1:Nn))/deltat;
    % FF=Mgm*aa(:)+Cgc*vv(:)+Kgk*xx(1:Nn);

    O_s2_st = obsv(Adkf2,H);
    
%% CHECKING THE OBSERVABILITY OF THE DAMAGES 
for z=1:Ndam
    Obsv_fd_s1(z) = Obsv_fd_s1(z) + norm(O_s1_fd(:, z))^2;
    
    Obsv_fd_s2(z) = Obsv_fd_s2(z) + norm(O_s2_fd(:, z))^2;
end

for k=1: 2*Nm

    Obsv_state_s1(k) = Obsv_state_s1(k) + norm(O_s1_st(:, k))^2;
    Obsv_state_s2(k) = Obsv_state_s2(k) + norm(O_s2_st(:, k))^2;
end
    
for k = 1: (2*Nm+Ndam)
    Obsv_ekf(k) = Obsv_ekf(k) + norm(O_ex(:,k))^2;

end
%     
%     

    
%    fuc(:,i)=inv(Hf)*(zf(:,i)-HBf);
    reba=1000;
    
   % fprintf(Fil1,'%f\n',fac');

      damfpredict(i,:)=fac;
   facc=fac;
%     
%     end
    P=(eye(dim)-K*H)*Pn;
    Pf=(eye(Ndam)-Kf*Hf)*Pnf;

     Pdkf2=(eye(dim)-Kx2*H)*Pndkf2;
    Pf2=(eye(Ndam)-Kdkf2*Hdkf2)*Pnf2;
   
    
xxp=xx;
vvp=vv;
aap=aa;

acp_2 = (xx2(Nm+1:2*Nm)-xxp2(Nm+1:2*Nm))/deltat;
xxp2 = xx2;


accp=acc;
   acp=ac;
   xp=x;
   vp=v;
   xrip=xri;
   vrip=vri;
   accpri=accri;

   
       
    
end
    
%% Vizualization 
% State plotting
figure(1);
sgtitle('Damage Estimation (sensor : z1,z3,z4,z5,Mid)');
for i=1:6
    subplot(2,3,i)
    plot(damf(i,:),'k-');
    hold on;
    plot(damfpredict(:,i),'b--'); hold on;
    plot(fd_est_dkf2(i,:),'m--');
    plot(fd_est_ex(i,:),'g--');
    xlabel('time');
    ylabel('Damage factor(fd)');
    hold off
end


% Damage factor plotting
figure(2);
for i=1:6
    subplot(2,3,i);
    plot(damf(i,:),'k-');
    hold on;
    plot(damfpredict(:,i),'b--'); hold on;
    plot(fd_est_dkf2(i,:),'m--');
    %plot(fd_est_ex(i,:),'g--');
    hold off;
end






% PLOTTING THE STATES  

figure(5);
  for i =1:Nm
        subplot(2,3,i);
        plot(linspace(0,deltat*N,N),xx_true_store(i,:),'k-'); hold on;
        plot(linspace(0,deltat*N,N),xx_store_dkf2(i,:),'m--');  hold on;
        plot(linspace(0,deltat*N,N),xx_store_dkf(i,:),'b--');
        xlabel('time');
        ylabel('displacement');
        title('disp(damage) vs time');

  end


%% Plotting the observation matrix 
figure(6)
bar(Obsv_fd_s1);%/max(Obsv_fd_s1));
set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6'});
title('Observability(Fd) over Time for scheme1');

figure(7)
bar(Obsv_fd_s2);%/max(Obsv_fd_s2));
set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6'});
title('Observability(Fd) over Time for scheme2 ');
figure(8)
bar(Obsv_state_s1);%/max(Obsv_state_s1));
set(gca, 'XTickLabel', {'m 1', 'm 2', 'm 3', 'm 4', 'm 5', 'm 6', 'm 7','m 8','m9','m10','m11','m12'});
title(' Observability Strength (modes) over Time for scheme1 ');

figure(9)
bar(Obsv_state_s2);%/max(Obsv_state_s2));
set(gca, 'XTickLabel', {'m 1', 'm 2', 'm 3', 'm 4', 'm 5', 'm 6', 'm 7','m 8','m9','m10','m11','m12'});
title(' Observability Strength (modes) over Time for scheme2 ');


figure(10)
bar(Obsv_ekf(2*Nm+1:2*Nm+Ndam));%/max(Obsv_ekf(2*Nm+1:2*Nm+Ndam)));
set(gca, 'XTickLabel', { 'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6'});
title(' Observability Strength (dams) over Time for EKF ');


P2_s2 = 1./Obsv_fd_s2;
P2_s1 = 1./Obsv_fd_s1;

P1_s1 = 1./Obsv_state_s1;
P1_s2 = 1./Obsv_state_s2;


P_ekf = 1./Obsv_ekf ; 

% printing 
for i =1:18
    P_ekf(i)
    %P2_s1(i)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         This the one with the 5 sensor 6 dam with the middle point of the zones   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% clear all
% syms tau
% p=1000;
% ro=2700; % ( in kg/m3)
% yo =70*10^9;%youngs modulus
% neu=0.3;
% h=0.001;%thickness
% a=0.65; %( m)
% b=3; % (m)
% divx=3;%input(' enter the no of divisions on the breadth ');
% % ex=length of element in x direction
% ex=a/(divx);
% divy=3;%input(' enter the no of divisions on the length ');
% % ey=length of element in y direction
% ey=b/(divy);
% 
% N=50000;
% deltat=0.00005;
% nObserv=6;
% 
% syms x y
% 
% %Fil1=fopen('facdat1.m','w');
% 
% 
% Nm=6;
% Nn=Nm
% %phi=[sin(pi*x/a)*sin(pi*y/b) sin(3*pi*x/a)*sin(pi*y/b) sinHafpi*x/a)*sin(3*pi*y/b) sin(3*pi*x/a)*sin(3*pi*y/b) sin(4*pi*x/a)*sin(4*pi*y/b)];
% phi=[sin(pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(2*pi*y/b) sin(2*pi*x/a)*sin(pi*y/b) sin(2*pi*x/a)*sin(2*pi*y/b) sin(3*pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(3*pi*y/b)];
% 
% dxxphi=diff(phi,x,2);
% dxyphi=diff(diff(phi,x),y);
% dyyphi=diff(phi,y,2);
% D=(yo*h^3)/(12*(1-neu^2));
% 
% Ndam=6;
% % damdox=[0 a/4;a/4 a/2; a/2 3*a/4;3*a/4 a;0 a/4;a/4 a/2];
% % damdoy=[0 b/4;0 b/4;0 b/4;0 b/4;b/4 b/2;b/4 b/2;b/4 b/2];
% 
% damdox=[0 a/3; a/3 2*a/3;2*a/3 a;0 a/3; a/3 2*a/3;2*a/3 a];
% damdoy=[0 b/2;0 b/2;0 b/2;b/2 b;b/2 b;b/2 b];
% 
% sensox=[a/6 a/2 5*a/6 a/6 a/2 5*a/6];
% sensoy=[b/4 b/4 b/4 3*b/4 3*b/4 3*b/4];
% 
% for i=1:Nm
%     for j=1:Nm
%         Mgm(i,j)=double(ro*h*int(int(phi(i)*phi(j),x,0,a),y,0,b));
%         Kgk(i,j)=double(D*int(int(dxxphi(i)*dxxphi(j)+2*dxyphi(i)*dxyphi(j)+dyyphi(i)*dyyphi(j),x,0,a),y,0,b));
%     end
%     Fgf(i,1)=double(int(int(p*phi(i),x,0,a),y,0,b));
% end
%    Cgc=0.0003*Mgm+0.0003*Kgk;     
% 
% for i=1:Ndam
%     for j=1:Nm
%         for k=1:Nm
%             Kd(i,j,k)=double(D*int(int(dxxphi(j)*dxxphi(k)+2*dxyphi(j)*dxyphi(k)+dyyphi(j)*dyyphi(k),x,damdox(i,1),damdox(i,2)),y,damdoy(i,1),damdoy(i,2)));
%         end
%     end
% end
% 
% 
% clear x
% gamma=1/2;
% beta=1/2;
% 
% x(:,1)=zeros(Nm,1);
% v(:,1)=zeros(Nm,1);
% acc(:,1)=zeros(Nm,1);
% ac(:,1)=zeros(Nm,1);
% Fg(:,1)=zeros(Nm,1);
% 
% Kgun=Kgk;
% Haf=[1 0 0 0];
% 
% damf=zeros(Ndam,N);
% % 
%   for i=2000:N
%       damf(1,i)=(i-2000)^.7*0.3/(N-2000)^.7;
%   end
% % 
% %   for i=12000:N
% %       damf(3,i)=(i-12000)*0.6/(N-12000);
% %   end
% % %  
%   for i=6000:N
%      damf(3,i)=(i-6000)*0.3/(N-6000);
%  end
% 
% for i=100:N
%     damf(4,i)=(i-100)*0.3/(N-100);
% end
% 
% for i=300:N
%     damf(5,i)=(i-300)*0.2/(N-300);
% end
% % 
% 
% 
% dim=2*Nm
% 
% % for i=1:nObserv
% %       Phi(i,:)=[phin1(Observ(i,1),Observ(i,2)) phin2(Observ(i,1),Observ(i,2)) phin3(Observ(i,1),Observ(i,2)) phin4(Observ(i,1),Observ(i,2))];
% % end
% 
% %Phi=[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
% 
% for i=1:nObserv
%     x1=sensox(i);
%     y1=sensoy(i);
%    % Phi(i,:)=[sin(pi*x1/a)*sin(pi*y1/b) sin(2*pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(2*pi*y1/b) sin(2*pi*x1/a)*sin(2*pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b) sin(3*pi*x1/a)*sin(3*pi*y1/b)  sin(pi*x1/a)*sin(4*pi*y1/b) sin(2*pi*x1/a)*sin(4*pi*y1/b) sin(3*pi*x1/a)*sin(4*pi*y1/b) sin(4*pi*x1/a)*sin(1*pi*y1/b) sin(4*pi*x1/a)*sin(2*pi*y1/b) sin(4*pi*x1/a)*sin(3*pi*y1/b)];
%     Phi(i,:)=[sin(pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(2*pi*y1/b) sin(2*pi*x1/a)*sin(pi*y1/b) sin(2*pi*x1/a)*sin(2*pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b)];
% 
% end
% %Phi=eye(Nm);
% 
% X=zeros(nObserv,Nm);
% 
% % H=[Phi X X;X Phi X;X X Phi];
% 
% H=[Phi X;X Phi];
% 
% Haf=Phi;
% 
% P=eye(dim);
% Pf=eye(Ndam)*10000;
% Q=eye(dim);
% %Q(2*Nn+1:dim,2*Nn+1:dim)=eye(dim-(2*Nn));
% dimf=1;
% 
% 
% Qf=diag([1 1 1 1 1 1]);
% 
% 
% % erR=z-[x;v];
% % erf=zf-acc;
% % R=cov(erR')*10000000;
% % Rf=cov(erf')*10000000;
% R=eye(2*nObserv)*10^(6);
% Rf=eye(nObserv)*10^(6);
% 
% xx=zeros(dim,1);
% xx(1,1)=0.0;
% 
% Fo(1)=0;
% 
% FF(:,1)=zeros(1,Nm);
% aa(:,1)=zeros(1,Nm);
% vv(:,1)=zeros(1,Nm);
% 
% fac=zeros(Ndam,1);
% facc=fac;
% 
% xxp=zeros(2*Nm,1)
% vvp=zeros(Nm,1);
% aap=zeros(Nm,1);
% accp=zeros(Nm,1);
%    acp=zeros(Nm,1);
%    xp=zeros(Nm,1);
%    vp=zeros(Nm,1);
%    xrip=zeros(Nm,1);
%    vrip=zeros(Nm,1);
%    accpri=zeros(Nm,1);
% 
% 
% for i=2:N
% 
%      Kgk=Kgun;
% %     
%     for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         Kgk=Kgk-damf(j,i)*tem;
%     end
% 
%     AA=Mgm+Cgc*deltat*gamma+Kgk*beta*deltat^2;
% 
% BB=Cgc*deltat*(1-gamma)+Kgk*(1-2*beta)*deltat^2/2;
% 
% CC=Cgc+deltat*Kgk;
% AAinv=inv(AA);
% 
%     t=(i-1)*deltat;
%     t
%     Fg(:)=Fgf*sin(2*pi/12*t);
%     acc(:)=AAinv*(Fg(:)-[BB]*accp(:)-[CC]*v(:)-Kgk*x(:));
%     accp=acc;
% 
%     for kk=1:Nm
%         ac(kk)=acc(kk)+0*(rand-.5)/50;
%         Accel(i,kk)=acc(kk);
%     end
%     accri=ac+(rand-.5)/70;
%     accri(1)=ac(1)+(rand-0.5)/5;
%     v(:)=vp+deltat*(1-gamma)*accpri+deltat*gamma*accri;
%     x(:)=xp+deltat*vp+(1-2*beta)/2*deltat^2*accpri+beta*deltat^2*accri;
%     XX(:,i)=x;
%     AC(:,i)=acc;
%     ACCRI(:,i)=accri;
%     vp=v;
%     xp=x;
% 
% 
%     for j=1:Nm
% 
%         vri(j)=vrip(j)+(acp(j)+ac(j))/2*deltat;
% 
%         xri(j)=xrip(j)+(vrip(j)+vri(j))/2*deltat;
%     end
% 
%     xrri=Phi*xri';
%     vrri=Phi*vri';
%    z(1:2*nObserv)=[xri(1:nObserv)'; vri(1:nObserv)'];
%   % z(1:2*nObserv,i)=[ v(1:nObserv,i) ; ac(1:nObserv,i)];
%    zf(1:nObserv)=Phi*ac;
% 
% 
% % end
% % 
% % 
% % 
% % xdiv=10;
% % ydiv=10;
% % for i=1:xdiv+1
% %     for j=1:ydiv+1
% %         xn=(i-1)*a/xdiv;
% %         yn=(j-1)*b/ydiv;
% %         phin1(i,j)=sin(pi*xn/a)*sin(pi*yn/b);
% %         phin2(i,j)=sin(3*pi*xn/a)*sin(pi*yn/b);
% %         phin3(i,j)=sin(pi*xn/a)*sin(3*pi*yn/b);
% %         phin4(i,j)=sin(3*pi*xn/a)*sin(3*pi*yn/b);
% %     end
% % end
% %         
% % xn=zeros(xdiv+1,ydiv+1,N);
% % vn=zeros(xdiv+1,ydiv+1,N);
% % an=zeros(xdiv+1,ydiv+1,N);
% % for i=1:(xdiv+1)
% %     for j=1:(ydiv+1)
% %       xn(i,j,:)=x(1,:)*phin1(i,j)+x(2,:)*phin2(i,j)+x(3,:)*phin3(i,j)+x(4,:)*phin4(i,j);
% %       vn(i,j,:)=v(1,:)*phin1(i,j)+v(2,:)*phin2(i,j)+v(3,:)*phin3(i,j)+v(4,:)*phin4(i,j);
% %       an(i,j,:)=ac(1,:)*phin1(i,j)+ac(2,:)*phin2(i,j)+ac(3,:)*phin3(i,j)+ac(4,:)*phin4(i,j);
% %     end
% % end
% % 
% % 
% %  
% % Observ=[6 6; 4 4];
% % 
% % % for j=1:nObserv
% % %   for i=1:N
% % %    z(j:nObserv:j+2*nObserv,i)=[xn(Observ(j,1),Observ(j,2),i)+0*(rand-0.5)*10^(-3)/5; vn(Observ(j,1),Observ(j,2),i)+0*(rand-0.5)*10^(-1)/3; an(Observ(j,1),Observ(j,2),i)];
% % %    zf(j,i)=an(Observ(j,1),Observ(j,2),i);
% % %   end
% % % end
% 
% 
% 
% 
% %for i=2:N
% 
%     iii=i
%     t=(i-1)*deltat;
% 
% 
% %     err=10000;
% %     
% %     
% %     while err>0.001
% %    
% %       
% 
%     disp=xxp(1:Nm)+xxp(Nm+1:2*Nm)*deltat;
%     vel=xxp(Nm+1:2*Nm)+acp(1:Nm)*deltat;
%     fac=facc;
%     Kdx=zeros(Nm,Ndam);
% 
%     for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         xi(:,j)=tem*disp;
%         Kdx(:,j)=xi(:,j);
%     end
% %        
%     F=Fgf*sin(2*pi/12*t);    
%     F1=F;
%     Mgminv=inv(Mgm);
% 
% 
%     Af=eye(Ndam);
%     Hf=Haf*Mgminv*Kdx;
%     Cgcx=Mgminv*(Kgun*disp+Cgc*vel);
%     HBf=Haf*(Mgminv*F-Cgcx);
% 
% 
% 
% 
%     Pnf=Af*Pf*transpose(Af)+Qf;
%     Kf=Pnf*transpose(Hf)*inv((Hf*Pnf*transpose(Hf))+Rf);
% 
%     Af=eye(1);
%     reba=100;
%     fac=facc+Kf*(zf'-Hf*fac(:)-HBf);
% 
% 
% 
%     Kgk1=Kgun;
%     for j=1:Ndam
%        for k=1:Nm
%            for m=1:Nm
%                tem(k,m)=Kd(j,k,m);
%            end
%        end
%        Kgk1=Kgk1-tem*fac(j);
%     end
% %     
% 
%     Ac=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgk1  -inv(Mgm)*Cgc];
%      Bc=[zeros(Nm,Nm); inv(Mgm)];
% % II=eye(Nn);
% % 
% % AA=Mgm+gamma*Cgc+beta*deltat^2*Kgk;
% % AAinv=inv(AA);
% % CC=Cgc+Kgk*deltat;
% % BB=deltat*Cgc*(1-gamma)+deltat^2/2*(1-2*beta)*Kgk;
% % 
% % A11=zeros(3*Nm,3*Nm);
% % 
% % A11(1:Nn,:)=[II-beta*deltat^2*AAinv*Kgk deltat*II-beta*deltat^2*AAinv*CC (1-2*beta)/2*deltat^2*II-beta*deltat^2*AAinv*BB];
% % 
% % A11(Nn+1:2*Nn,:)=[-gamma*deltat*AAinv*Kgk II-gamma*deltat*AAinv*CC (1-gamma)*deltat*II-gamma*deltat*AAinv*BB];
% % 
% % A11(2*Nn+1:3*Nn,:)=AAinv*[-Kgk -CC -BB];
% % A=A11;
% % 
% % B=[eye(Nm,Nm); gamma*deltat*eye(Nm,Nm); beta*deltat^2*eye(Nm,Nm)];
% 
% % A=exp(Ac*deltat);
% % 
% % B=double(int(exp(Ac*tau),tau,0,deltat))*Bc;
% 
% A=deltat*Ac+eye(2*Nm);
% B=deltat*Bc;
% 
% A=zeros(2*Nm,2*Nm);
% BB=zeros(2*Nm,2*Nm);
% B=zeros(2*Nm,Nm);
% for jj=0:10
%     A=A+(Ac*deltat)^jj/factorial(jj);
%     BB=BB+(Ac*deltat)^jj/factorial(jj+1);
% end
% 
% B=BB*Bc*deltat;
% 
%      %A=exp(Ac*deltat);
%      %B=[A-eye(2*Nn)]*inv(Ac)*Bc;
% %      A=inv(eye(2*Nm)-Ac);
% %      B=inv(eye(2*Nm)-Ac)*Bc;
% 
% 
% 
%    % 0*B*Fg(:,i-1)
%     xx=A*xxp+B*F;
%     Pn=A*P*transpose(A)+Q;    
%     K=Pn*transpose(H)*inv((H*Pn*transpose(H))+R);
%     xx=xx+K*(z'-H*xx);    
%     vv=(xx(1:Nn)-xxp(1:Nn))/deltat;
%     aa=(vv(1:Nn)-vvp(1:Nn))/deltat;
%     FF=Mgm*aa(:)+Cgc*vv(:)+Kgk*xx(1:Nn);
%     %KK(:,i)=K;
% 
% 
% 
%     %Fo(i)=Fo(i-1)+deltat/2*(xx(5,i-1)+xx(5,i));
% 
% 
% 
% 
% %     
% %     
% 
% 
% %    fuc(:,i)=inv(Hf)*(zf(:,i)-HBf);
%     reba=1000;
% 
%    % fprintf(Fil1,'%f\n',fac');
% 
%       damfpredict(i,:)=fac;
%    facc=fac;
% %     
% %     end
%     P=(eye(dim)-K*H)*Pn;
%     Pf=(eye(Ndam)-Kf*Hf)*Pnf;
% 
% 
% xxp=xx;
% vvp=vv;
% aap=aa;
% accp=acc;
%    acp=ac;
%    xp=x;
%    vp=v;
%    xrip=xri;
%    vrip=vri;
%    accpri=accri;
% 
% 
% 
% 
% end
% 
% 
% for i=1:6
%     subplot(2,3,i)
%     plot(damf(i,:))
%     hold on
%     plot(damfpredict(:,i))
%     hold off
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Different set of code    
%
%---------------------------------------------------------------------------------------------------------------------------------
% %% Structural Health Monitoring: Stabilized Modal EKF
% clearvars; clc; close all;
% syms x y
% 
% % --- 1. Physical Parameters ---
% a = 0.6; b = 0.4; h = 0.001; 
% ro = 2700; yo = 70e9; neu = 0.3;
% p_load = 1000; 
% D = (yo*h^3)/(12*(1-neu^2));
% 
% % --- 2. Simulation & Modal Parameters ---
% Nm = 6; Ndam = 6; nObserv = 6;
% deltat = 0.00005; % Time step
% N = 15000; % Total steps
% w_drive = 137.2892; 
% 
% % Define Mode Shapes
% phi = [sin(pi*x/a)*sin(pi*y/b), sin(pi*x/a)*sin(2*pi*y/b), ...
%        sin(2*pi*x/a)*sin(pi*y/b), sin(2*pi*x/a)*sin(2*pi*y/b), ...
%        sin(pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(1*pi*y/b)];
% 
% dxxphi = diff(phi,x,2); dxyphi = diff(diff(phi,x),y); dyyphi = diff(phi,y,2);
% 
% % --- 3. System Matrix Assembly ---
% fprintf('Assembling Matrices...\n');
% Mgm = zeros(Nm); Kgk = zeros(Nm); Fgf = zeros(Nm,1);
% for i=1:Nm
%     for j=1:Nm
%         Mgm(i,j) = double(ro*h*int(int(phi(i)*phi(j),x,0,a),y,0,b));
%         Kgk(i,j) = double(D*int(int(dxxphi(i)*dxxphi(j)+2*dxyphi(i)*dxyphi(j)+dyyphi(i)*dyyphi(j),x,0,a),y,0,b));
%     end
%     Fgf(i,1) = double(int(int(p_load*phi(i),x,0,a),y,0,b));
% end
% Cgc = 0.0003*Mgm + 0.0003*Kgk;
% 
% % --- 4. Mass Normalization (CRITICAL FOR STABILITY) ---
% % We transform the system so Mass Matrix = Identity
% [V, D_eig] = eig(Kgk, Mgm);
% T = V; % Transformation matrix
% Mgm_n = T' * Mgm * T; 
% Kgk_n = T' * Kgk * T;
% Cgc_n = T' * Cgc * T;
% Fgf_n = T' * Fgf;
% 
% % Re-scale Kd for normalized space
% Kd = zeros(Ndam, Nm, Nm);
% damdox = [0 a/3; a/3 2*a/3; 2*a/3 a; 0 a/3; a/3 2*a/3; 2*a/3 a];
% damdoy = [0 b/2; 0 b/2; 0 b/2; b/2 b; b/2 b; b/2 b];
% for i=1:Ndam
%     for j=1:Nm
%         for k=1:Nm
%             val = double(D*int(int(dxxphi(j)*dxxphi(k)+2*dxyphi(j)*dxyphi(k)+dyyphi(j)*dyyphi(k),x,damdox(i,1),damdox(i,2)),y,damdoy(i,1),damdoy(i,2)));
%             Kd(i,j,k) = val;
%         end
%     end
%     % Transform each Kd zone to normalized coordinates
%     temp_Kd = squeeze(Kd(i,:,:));
%     Kd(i,:,:) = T' * temp_Kd * T;
% end
% 
% % Sensor locations and Phi matrix
% sensox = [a/6 a/2 5*a/6 a/6 a/2 5*a/6];
% sensoy = [b/4 b/4 b/4 3*b/4 3*b/4 3*b/4];
% Phi_orig = zeros(nObserv, Nm);
% for i=1:nObserv
%     for j=1:Nm
%         Phi_orig(i,j) = double(subs(phi(j), {x,y}, {sensox(i), sensoy(i)}));
%     end
% end
% Phi = Phi_orig * T; % Normalized Phi
% 
% % --- 5. EKF Initialization ---
% nx = 2*Nm + Ndam;
% X_hat = zeros(nx, 1);
% P = diag([ones(Nm,1)*1e-6; ones(Nm,1)*1e-4; ones(Ndam,1)*1e-2]); 
% Q = diag([ones(Nm,1)*1e-9; ones(Nm,1)*1e-7; ones(Ndam,1)*1e-5]);
% R = 1e-1 * eye(nObserv); % Increased noise floor for stability
% 
% % Simulation Storage
% true_theta = zeros(Ndam, N);
% for i=3000:N, true_theta(1,i) = 0.25; end % Step damage
% theta_est = zeros(Ndam, N);
% q_t = zeros(Nm,1); qd_t = zeros(Nm,1);
% 
% % --- 6. Main Loop ---
% for k = 2:N
%     t = (k-1)*deltat;
%     F_in = Fgf_n * sin(w_drive * t);
% 
%     % --- Plant Simulation (Normalized) ---
%     K_now = Kgk_n;
%     for d = 1:Ndam
%         K_now = K_now - true_theta(d,k) * squeeze(Kd(d,:,:));
%     end
%     qdd_t = F_in - Cgc_n*qd_t - K_now*q_t; % M is now Identity
%     qd_t = qd_t + qdd_t * deltat;
%     q_t = q_t + qd_t * deltat;
%     zf = Phi * qdd_t + 0.01*randn(nObserv,1); % Simulated acceleration
% 
%     % --- EKF Prediction ---
%     q_e = X_hat(1:Nm); qd_e = X_hat(Nm+1:2*Nm); th_e = X_hat(2*Nm+1:end);
% 
%     K_e = Kgk_n; S_mat = zeros(Nm, Ndam);
%     for d = 1:Ndam
%         Ki = squeeze(Kd(d,:,:));
%         K_e = K_e - th_e(d) * Ki;
%         S_mat(:,d) = Ki * q_e;
%     end
% 
%     % State Space Jacobian (A)
%     Ac = [zeros(Nm), eye(Nm), zeros(Nm, Ndam);
%           -K_e, -Cgc_n, S_mat;
%           zeros(Ndam, nx)];
% 
%     % Discrete Transition (Exponential Map)
%     Ax = expm(Ac * deltat); 
% 
%     X_pred = Ax * X_hat;
%     P_pred = Ax * P * Ax' + Q;
% 
%     % Observation Jacobian (H) - Acceleration
%     H = [-Phi*K_e, -Phi*Cgc_n, Phi*S_mat];
% 
%     % --- Innovation & Correction ---
%     h_x = Phi * (F_in - Cgc_n*qd_e - K_e*q_e);
%     Inn = zf - h_x;
%     S_cov = H * P_pred * H' + R;
% 
%     % NAN GUARD: Check Condition Number
%     if rcond(S_cov) < 1e-15 || any(isnan(S_cov), 'all')
%         fprintf('Warning: S_cov singular at k=%d. Skipping update.\n', k);
%         X_hat = X_pred; P = P_pred;
%     else
%         K_gain = (P_pred * H') / S_cov;
%         % Limit the gain to prevent explosion
%         K_gain = max(min(K_gain, 100), -100); 
%         X_hat = X_pred + K_gain * Inn;
%         P = (eye(nx) - K_gain * H) * P_pred;
%     end
% 
%     % Final Safety
%     X_hat(isnan(X_hat)) = 0; 
%     theta_est(:,k) = X_hat(2*Nm+1:end);
% end
% 
% % --- 7. Plotting ---
% figure;
% for i=1:Ndam
%     subplot(3,2,i);
%     plot(true_theta(i,:), 'k'); hold on;
%     plot(theta_est(i,:), 'r--');
%     title(['Zone ', num2str(i)]); ylim([-0.1 0.5]);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4s 6 dam

%% Done modal domain with observability 


clear all

p=1000;
ro=2700; % ( in kg/m3)
yo =70*10^9;%youngs modulus
neu=0.3;
h=0.001;%thickness
a=0.6; %( m)
b=0.4; % (m)
divx = 3;  %input(' enter the no of divisions on the breadth ');
% ex=length of element in x direction
ex=a/(divx);
divy=3;  %input(' enter the no of divisions on the length ');
% ey=length of element in y direction
ey=b/(divy);
%% Natural frequencies for the above 
%% NB: when frequencies are around 8,9 w1 we find scheme 2 is better than the other model
w_1 = 100*pi; %1*137.2892;
w_6 =  897.6603; 

N = 5000;
deltat = 0.00005;
nObserv= 4 ;       % No of observable zones (** sensor is acting) 

syms x y




Nm=6;
Nn=Nm;
%phi=[sin(pi*x/a)*sin(pi*y/b) sin(3*pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(3*pi*y/b) sin(3*pi*x/a)*sin(3*pi*y/b) sin(4*pi*x/a)*sin(4*pi*y/b)];
phi=[sin(pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(2*pi*y/b) sin(2*pi*x/a)*sin(pi*y/b) sin(2*pi*x/a)*sin(2*pi*y/b) sin(1*pi*x/a)*sin(3*pi*y/b) sin(3*pi*x/a)*sin(1*pi*y/b)];
%phi=[sin(pi*x/a)*sin(pi*y/b)  sin(3*pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(3*pi*y/b)  sin(3*pi*x/a)*sin(3*pi*y/b)  sin(pi*x/a)*sin(5*pi*y/b)  sin(5*pi*x/a)*sin(1*pi*y/b)];

dxxphi=diff(phi,x,2);
dxyphi=diff(diff(phi,x),y);
dyyphi=diff(phi,y,2);
D=(yo*h^3)/(12*(1-neu^2));

Ndam=6;  % Total no. of damage zones ( each having own fd) not the sensors 
% damdox=[0 a/4;a/4 a/2; a/2 3*a/4;3*a/4 a;0 a/4;a/4 a/2];
% damdoy=[0 b/4;0 b/4;0 b/4;0 b/4;b/4 b/2;b/4 b/2;b/4 b/2];

damdox=[0 a/3; a/3 2*a/3;2*a/3 a;0 a/3; a/3 2*a/3;2*a/3 a];
damdoy=[0 b/2;0 b/2;0 b/2;b/2 b;b/2 b;b/2 b];

% Most obsv Sensor location for 6 sensor placement 
% sensox = [0.2949  0.1627  0.4271  0.1525  0.4475  0.3051];
% sensoy = [ 0.1949 0.1026  0.1026  0.2769  0.2769  0.3179];


% sensox = [0.2949  0.1627  0.4271  0.1525  0.4475 ];
% sensoy = [ 0.1949 0.1026  0.1026  0.2769  0.2769 ];

% Middle term is added as the 
 sensox=[a/6  a/2  5*a/6  a/6   5*a/6   a/2 ];
 sensoy=[b/4  b/4  b/4   3*b/4  3*b/4   b/2 ];


 % % Effective for 5 sensor location 
 % sensox=[  a/6   5*a/6   1*a/6   5*a/6    a/2   ];
 % sensoy=[   b/4   b/4    3*b/4    3*b/4    b/2 ];

 % % Effective for 3 sensor location 
 % sensox=[  a/6    5*a/6    a/2   ];
 % sensoy=[   b/4    3*b/4    b/2 ];

 %  % Effective for 2 sensor location 
 % sensox=[  a/6        a/2   ];
 % sensoy=[   b/4       b/2 ];



for i=1:Nm
    for j=1:Nm
        Mgm(i,j)=double(ro*h*int(int(phi(i)*phi(j),x,0,a),y,0,b));
        Kgk(i,j)=double(D*int(int(dxxphi(i)*dxxphi(j)+2*dxyphi(i)*dxyphi(j)+dyyphi(i)*dyyphi(j),x,0,a),y,0,b));
    end
    Fgf(i,1)=double(int(int(p*phi(i),x,0,a),y,0,b));
end
   Cgc=0.0003*Mgm+0.0003*Kgk;     

for i=1:Ndam
    for j=1:Nm
        for k=1:Nm
            Kd(i,j,k)=double(D*int(int(dxxphi(j)*dxxphi(k)+2*dxyphi(j)*dxyphi(k)+dyyphi(j)*dyyphi(k),x,damdox(i,1),damdox(i,2)),y,damdoy(i,1),damdoy(i,2)));
        end
    end
end
    

clear x
gamma=1/2;
beta=1/2;

x(:,1)=zeros(Nm,1);
v(:,1)=zeros(Nm,1);
acc(:,1)=zeros(Nm,1);
ac(:,1)=zeros(Nm,1);
Fg(:,1)=zeros(Nm,1);

Kgun=Kgk;  % kGK Is the undamaged one 
Haf=[1 0 0 0];

damf=zeros(Ndam,N);
% 
  for i=200:N
      damf(1,i)=(i-200)^.7*0.3/(N-200)^.7;
  end
% % 
  for i=1200:N
      damf(2,i)=(i-1200)*0.6/(N-1200);
  end

%   for i=600:N
%      damf(3,i)=(i-600)*0.12/(N-600);
%  end
% 
% for i=10:N
%     damf(4,i)=((i-10)^0.99)*0.3/(N-10);
% end

for i=30:N
    damf(5,i)=(i-30)*0.2/(N-30);
end


for i=300:N
    damf(6,i)=(i-300)*0.2/(N-300);
end

% % NON LINEAR DAMAGES 
% for i=2000:N
%       damf(1,i)=(i-2000)^.7*0.3/(N-2000)^.7;
%   end
% % % 
%   for i=12000:N
%       damf(2,i)=(i-12000)*0.8/(N-12000);
%   end
% %  
%   for i=600:N
%      damf(3,i)=((i-600)^0.8)*0.12/(N-600);
%  end
% 
% for i=100:N
%     damf(4,i)=(i-100)*0.3/(N-100);
% end
% 
% for i=300:N
%     damf(5,i)=(i-300)*0.2/(N-300);
% end
% 
% 
% for i=300:N
%     damf(6,i)=(i-300)*0.2/(N-300);
% end
% 


dim=2*Nm;

% for i=1:nObserv
%       Phi(i,:)=[phin1(Observ(i,1),Observ(i,2)) phin2(Observ(i,1),Observ(i,2)) phin3(Observ(i,1),Observ(i,2)) phin4(Observ(i,1),Observ(i,2))];
% end

%Phi=[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];

for i=1:nObserv
    x1=sensox(i);
    y1=sensoy(i);
   % Phi(i,:)=[sin(pi*x1/a)*sin(pi*y1/b) sin(2*pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(2*pi*y1/b) sin(2*pi*x1/a)*sin(2*pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b) sin(3*pi*x1/a)*sin(3*pi*y1/b)  sin(pi*x1/a)*sin(4*pi*y1/b) sin(2*pi*x1/a)*sin(4*pi*y1/b) sin(3*pi*x1/a)*sin(4*pi*y1/b) sin(4*pi*x1/a)*sin(1*pi*y1/b) sin(4*pi*x1/a)*sin(2*pi*y1/b) sin(4*pi*x1/a)*sin(3*pi*y1/b)];
     %Phi(i,:) = [sin(pi*x1/a)*sin(pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b)  sin(3*pi*x1/a)*sin(3*pi*y1/b)  sin(pi*x1/a)*sin(5*pi*y1/b)  sin(5*pi*x1/a)*sin(1*pi*y1/b)];
     Phi(i,:) = [ sin(pi*x1/a)*sin(pi*y1/b) sin(1*pi*x1/a)*sin(2*pi*y1/b) sin(2*pi*x1/a)*sin(1*pi*y1/b) sin(2*pi*x1/a)*sin(2*pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b)];
end
maxVal = max(Phi, [], 'all');
Phi = Phi./maxVal;


X=zeros(nObserv,Nm);


% H=[Phi X X;X Phi X;X X Phi];

H=[Phi  X ; X  Phi];

Haf=Phi;
% 
% 
 P=[ones(Nm,1)*10^(-2);ones(Nm,1)*10^(-2)];   % 10^-10  for in 4 sensors
 P=diag(P);
 Pdkf2=P;
%P = 1e-5*((3.8302e-61)^-1)*blkdiag(1.3013e-65,2.6286e-69, 3.4733e-66, 2.8234e-70, 1.6619e-73,6.9362e-69 , 3.2922e-61 , 6.7425e-66, 3.8302e-61, 3.1583e-65, 5.8089e-70, 1.1842e-63);   % for 6 sensor and usual mode shapes
%Pdkf2= 1e-5*((3.8302e-61)^-1)*blkdiag( 1.3047e-65,2.6293e-69, 3.4802e-66, 2.8268e-70, 1.6648e-73,6.9442e-69 , 3.2971e-61 , 6.7527e-66, 3.8364e-61, 3.1618e-65, 5.8148e-70, 1.1854e-63); % for 6 sensor and usual mode shapes

% % % --- 5 sensor based 6 Ndam  100 pi freq ---------
Pf = 1e-8*blkdiag( 0.1073  , 0.0543 , 0.1287,0.1229,  0.0521 , 0.1145);% for 6 sensor scheme 1;
Pf2 = 1e-8*blkdiag(0.1077 , 0.0544 ,0.1292, 0.1234 ,  0.0522,  0.1149); % for 6 sensor scheme 2


 
% Pf = 1e-8*blkdiag(0.1217,0.0553,0.1344,0.1371, 0.0976, 0.4554); % for 4 sensor 
% Pf2 = 1e-8*blkdiag(0.1221,0.0554,0.1347,0.1376, 0.0982 ,  0.4585);


 Q=eye(dim)/10^-10;  %  % scheme 1 -> eye(dim)/10^-10 for 6 sensor (unchange modes)
 Qdkf2=eye(dim)/10^-10; % scheme 2 -> eye(dim)/10^-8 for 6 sensor (unchange modes)
%Q = 1e5*P; % for 6 sensor unchange modes
%Qdkf2=1e5*Pdkf2;
%Q(2*Nn+1:dim,2*Nn+1:dim)=eye(dim-(2*Nn));
dimf=1;

%Pekf = 1e-10*eye(2*Nm + Ndam);  %working well in 6s (usual modes) 
%Pekf = 1e6*blkdiag(3.3898e-15,3.8276e-17, 2.4722e-16, 1.3198e-17, 1.8844e-18,2.3570e-17, 3.4714e-08 ,3.8733e-10, 2.5383e-09,  1.3562e-10, 1.8996e-11, 2.4181e-10, 1.9583e-13, 1.3443e-13 ,  2.0059e-13, 1.9658e-13, 1.3482e-13, 2.0093e-13  ); % gives the modes 1e1*multiple
%Qekf = 1e-8*eye(2*Nm + Ndam);
Q_ekf_fd = 1e-2*blkdiag(0.0249,0.0407,0.0346, 0.0283, 0.266 , 0.8991 );      %    0.2663 , 0.2135 , 0.3137, 0.1538 ,0.4348, 0.8606);
P_state_ekf = 1e0*blkdiag(  2.3679e-15 , 4.8768e-17 , 2.1346e-16 , 5.6883e-18 , 1.2955e-18, 9.1274e-18, 2.5246e-08 , 5.1644e-10, 2.2563e-09, 6.0273e-11, 1.3815e-11, 9.6723e-11 );
%Pekf = blkdiag(1e8*P_state_ekf,Q_ekf_fd);

Pekf = 1e2*blkdiag( 8.1289e-16 ,3.7103e-16,1.6217e-15,4.3337e-17,1.2132e-17,7.7268e-17,  2.3118e-07, 3.8691e-09, 1.6926e-08 ,4.5216e-10, 1.2650e-10 , 8.0516e-10,9.6168e-09 ,4.4077e-09, 9.6275e-09,9.1041e-09, 4.3114e-09, 9.2748e-09);
Qekf = 1e-6*Pekf; %blkdiag(3.3898e-15,3.8276e-17, 2.4722e-16, 1.3198e-17, 1.8844e-18,2.3570e-17, 3.4714e-08 ,3.8733e-10, 2.5383e-09,  1.3562e-10, 1.8996e-11, 2.4181e-10, 1.9583e-13, 1.3443e-13 ,  2.0059e-13, 1.9658e-13, 1.3482e-13, 2.0093e-13  );


%Qf=    ;
%Qf = 1e-11*blkdiag( 0.2438  ,0.1030 ,0.2463  ,0.2365   , 0.1000, 0.2328);  % working better in 1e-12
Qf =  1e0*Pf; %1e-12*blkdiag( 0.4793  , 0.3464 ,0.5248, 0.4994 ,0.3564, 0.5304);
Qf2 = 1e0*Pf2; %    % 1e1*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);


R=1*eye(2*nObserv)*10^(2);  % for state ( disp & vel) of Scheme - 1 for 6 sen ->1*eye(2*nObserv)*10^(-4) works better ( no noise)  
Rf=1*eye(nObserv)*10^(0.2);
Rf2=1*eye(nObserv)*10^(0.2);
%Rf = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);

Rdkf2=1*eye(2*nObserv)*10^(2); % for state ( disp & vel) of Scheme - 2
%Rf2 = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);

Rekf = 10^(0.1)*eye(nObserv);
% -------------------------------------------  5s --------------------------- dsada-------------------------j







%% nO OF THE sensor reduced Increase the R value

xx=zeros(dim,1);
xx(1,1)=0.0;

Fo(1)=0;

FF(:,1)=zeros(1,Nm);
aa(:,1)=zeros(1,Nm);
vv(:,1)=zeros(1,Nm);

fac=zeros(Ndam,1);
facc=fac;
fd_ex = zeros(Ndam,1);
fd_est_ex = zeros(Ndam,N);
fd_dkf2 = zeros(Ndam,1);
fd_est_dkf2 = zeros(Ndam,N);

xxp=zeros(2*Nm,1);
xxp2=zeros(2*Nm,1);
xxp_exfd=zeros(2*Nm+Ndam,1);  % for the extended kalman filter 
vvp=zeros(Nm,1);
aap=zeros(Nm,1);
accp=zeros(Nm,1);
   acp=zeros(Nm,1);
   acp_2 = zeros(Nm,1);
   xp=zeros(Nm,1);
   vp=zeros(Nm,1);
   xrip=zeros(Nm,1);
   vrip=zeros(Nm,1);
   accpri=zeros(Nm,1);
xx_true_store = zeros(Nm,N);   
xx_store_dkf = zeros(Nm,N);
xx_store_dkf2 = zeros(Nm,N);
xx_store_ekf = zeros(Nm,N);

%% Creating some obs matrix for the different schemes and the states ( sc1 & sc2 & extened) 
Obsv_state_s1 = zeros(2*Nm,1);
Obsv_fd_s2 = zeros(Ndam,1);
Obsv_fd_s1 = zeros(Ndam,1);
Obsv_state_s2 = zeros(2*Nm,1);
Obsv_ekf = zeros(2*Nm + Ndam, 1);

for i=2:N
    
     Kgk=Kgun;
%     
    for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgk=Kgk-damf(j,i)*tem;
    end
   
    AA=Mgm+Cgc*deltat*gamma+Kgk*beta*deltat^2;

BB=Cgc*deltat*(1-gamma)+Kgk*(1-2*beta)*deltat^2/2;

CC=Cgc+deltat*Kgk;
AAinv=inv(AA);
    
    t=(i-1)*deltat;
    t;
    Fg(:)=Fgf*sin(w_1*t);
    acc(:)=AAinv*(Fg(:)-[BB]*accp(:)-[CC]*v(:)-Kgk*x(:));
    accp=acc;
    
    for kk=1:Nm
        ac(kk)=acc(kk)+ 0*(rand-.5)/50;
        Accel(i,kk)=acc(kk);
    end
    accri=ac+ 0*(rand-.5)/70;
    accri(1)=ac(1)+ 0*(rand-0.5)/5;
    v(:)=vp+deltat*(1-gamma)*accpri+deltat*gamma*accri;
    x(:)=xp+deltat*vp+(1-2*beta)/2*deltat^2*accpri+beta*deltat^2*accri;
    XX(:,i)=x;
    AC(:,i)=acc;
    ACCRI(:,i)=accri;
    vp=v;
    xp=x;
    
    
    for j=1:Nm
        
        vri(j)=vrip(j)+(acp(j)+ac(j))/2*deltat;
        
        xri(j)=xrip(j)+(vrip(j)+vri(j))/2*deltat;
    end

    xx_true_store(:,i) = xri';
    
    xrri=Phi*(xri+xri*0.002*rand)';
    vrri=Phi*(vri+vri*0.002*rand)';
    z = [xrri; vrri]';
   %z(1:2*nObserv)=[xri(1:nObserv)'; vri(1:nObserv)'];
  % z(1:2*nObserv,i)=[ v(1:nObserv,i) ; ac(1:nObserv,i)];
   zf(1:nObserv)=Phi*(ac+0.002*rand);

   





%for i=2:N
    
    iii=i;
    t=(i-1)*deltat;
    
   

    disp=xxp(1:Nm)+xxp(Nm+1:2*Nm)*deltat;
    vel=xxp(Nm+1:2*Nm)+acp(1:Nm)*deltat;
    fac=facc;
    Kdx=zeros(Nm,Ndam);
  
    for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        xi(:,j)=tem*disp;
        Kdx(:,j)=xi(:,j);
    end
%       F=Fgf*sin(2*pi*t/6);
    F=Fgf*sin(w_1*t);    
    F1=F;
    Mgminv=inv(Mgm);
    
    
    Af=eye(Ndam);
    Hf=Haf*Mgminv*Kdx;
    Cgcx=Mgminv*(Kgun*disp+Cgc*vel);
    HBf=Haf*(Mgminv*F-Cgcx);

    O_s1_fd = obsv(Af,Hf);
    
    % for Scheme 2 

     Kgks2 = Kgun;
     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgks2=Kgks2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
     end

      Kdx2=zeros(Nm,Ndam);
  
    for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        xi(:,j)=tem*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2 );
        Kdx2(:,j)=xi(:,j);
    end
    Hdkf2 = Phi*(inv(Mgm + deltat*0.5*Cgc + deltat^2*0.25*(Kgks2))*(Kdx2)) ;
    a_dkf2 = Hdkf2*fd_dkf2 + Phi*(inv(Mgm + deltat*0.5*Cgc + deltat^2*0.25*(Kgks2))*(F1- Kgun*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2) - Cgc*(xxp2(Nm+1:2*Nm) + deltat^2*0.5*acp_2) ));
    
   O_s2_fd = obsv(Af,Hdkf2);
    % for scheme 1 
    Pnf=Af*Pf*transpose(Af)+Qf;
    Kf=Pnf*transpose(Hf)*inv((Hf*Pnf*transpose(Hf))+Rf);

    %Af=eye(1);
    reba=100;
    fac=facc+Kf*(zf'-Hf*fac(:)-HBf);

     %for Scheme 2
     Pnf2=Af*Pf2*transpose(Af)+Qf2;
    Kdkf2=Pnf2*transpose(Hdkf2)*inv((Hdkf2*Pnf2*transpose(Hdkf2))+Rf2);

    %Af=eye(1);
    reba=100;
    fd_dkf2=fd_dkf2+Kdkf2*(zf'-a_dkf2);






%% 2nd Loop for DEKF 

    Kgk1 = Kgun;
     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgk1=Kgk1-fac(j)*tem;  % here must be fac or facc damf will not be there
     end


    Ac=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgk1  -inv(Mgm)*Cgc];
     Bc=[zeros(Nm,Nm); inv(Mgm)];

 % FOr Scheme 2 

    Kgdkf2 = Kgun;
     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgdkf2=Kgdkf2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
     end


    Adkf2=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgdkf2  -inv(Mgm)*Cgc];
     Bdkf2=[zeros(Nm,Nm); inv(Mgm)];




% -------- Extended kalman filter matrix --------
     S_ekf=zeros(Nm,Ndam);
    for j=1:Ndam
       for k=1:Nm
           for m=1:Nm
               tem(k,m)=Kd(j,k,m);
           end
       end
       S_ekf(:,j) = tem*xxp_exfd(1:Nm);
    end
 
    Aex = [zeros(Nm,Nm) , eye(Nm,Nm) , zeros(Nm,Ndam);
            -inv(Mgm)*Kgk1 , -inv(Mgm)*Cgc , inv(Mgm)*S_ekf ;
            zeros(Ndam,Nm), zeros(Ndam,Nm), eye(Ndam,Ndam)];

    Bex = [zeros(Nm,Nm); inv(Mgm); zeros(Ndam,Nm)];


% Discretization 
% for scheme 1
A=deltat*Ac+eye(2*Nm);
B=deltat*Bc;

% How this is defined
A=zeros(2*Nm,2*Nm);
BB=zeros(2*Nm,2*Nm);
B=zeros(2*Nm,Nm);
for jj=0:10
    A=A+(Ac*deltat)^jj/factorial(jj);
    BB=BB+(Ac*deltat)^jj/factorial(jj+1);
end

B=BB*Bc*deltat;



     %A=exp(Ac*deltat);
     %B=[A-eye(2*Nn)]*inv(Ac)*Bc;
%      A=inv(eye(2*Nm)-Ac);
%      B=inv(eye(2*Nm)-Ac)*Bc;

% for scheme 2 
Ax2=deltat*Adkf2+eye(2*Nm);
Bx2=deltat*Bdkf2;

% How this is defined
Ax2=zeros(2*Nm,2*Nm);
Bx2_=zeros(2*Nm,2*Nm);
Bx2=zeros(2*Nm,Nm);
for jj=0:10
    Ax2=Ax2+(Adkf2*deltat)^jj/factorial(jj);
    Bx2_=Bx2_+(Adkf2*deltat)^jj/factorial(jj+1);
end

Bx2=Bx2_*Bdkf2*deltat;


% ----- For extended kalman filter 
    Aex =  deltat*Aex + eye(2*Nm + Ndam);
    Bex=deltat*Bex;

    Ax=zeros(2*Nm+Ndam,2*Nm+Ndam);
    BBx=zeros(2*Nm+Ndam,2*Nm+Ndam);
    Bx=zeros(2*Nm+Ndam,Nm);
    for jj=0:10
        Ax=Ax+(Aex*deltat)^jj/factorial(jj);
        BBx=BBx+(Aex*deltat)^jj/factorial(jj+1);
    end

    Bx=BBx*Bex*deltat;

    xx_exfd = Ax*xxp_exfd + Bx*F;
    Pekfn = Ax*Pekf*transpose(Ax)+ Qekf; 
     Kgk2 = Kgun;

     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgk2=Kgk2-xx_exfd(2*Nm+j)*tem;  % here must be fac or facc damf will not be there
    end
      S_ekf=zeros(Nm,Ndam);
    for j=1:Ndam
       for k=1:Nm
           for m=1:Nm
               tem(k,m)=Kd(j,k,m);
           end
       end
       S_ekf(:,j) = tem*xx_exfd(1:Nm);
    end
 


    Hekf = [-Phi*inv(Mgm)*Kgk2 , -Phi*inv(Mgm)*Cgc , Phi*inv(Mgm)*S_ekf ];
    z_diff_ekf = zf' - Phi*( Mgm\(F1 - Kgk2*xx_exfd(1:Nm) - Cgc*xx_exfd(Nm+1 :2*Nm) )) ;
    Kekf = Pekfn*transpose(Hekf)*inv((Hekf*Pekfn*transpose(Hekf)) + Rekf);
    xxp_exfd = xx_exfd + Kekf*z_diff_ekf;

    % Storing the Values of the updated ekf 
    xx_store_ekf(:,i) = xxp_exfd(1:Nm); %(  first Nm 'u' then Nm 'vel')
    fd_est_ex(:,i) = xxp_exfd(2*Nm +1: 2*Nm + Ndam);

    O_ex = obsv(Ax,Hekf);

 
   % 0*B*Fg(:,i-1)
   % For Scheme 1 
    xx=A*xxp+B*F;
    Pn=A*P*transpose(A)+Q;    
    K=Pn*transpose(H)*inv((H*Pn*transpose(H))+R);
    xx=xx+K*(z'-H*xx);    
    xx_store_dkf(:,i)=xx(1:Nm);
    vv=(xx(1:Nn)-xxp(1:Nn))/deltat;
    aa=(vv(1:Nn)-vvp(1:Nn))/deltat;
    FF=Mgm*aa(:)+Cgc*vv(:)+Kgk*xx(1:Nn);
    %KK(:,i)=K;
     O_s1_st = obsv(Ac,H);
    
    
    %Fo(i)=Fo(i-1)+deltat/2*(xx(5,i-1)+xx(5,i));
    
    % For Scheme 2 
    xx2=Ax2*xxp2+Bx2*F;
    Pndkf2=Ax2*Pdkf2*transpose(Ax2)+Qdkf2;    
    Kx2=Pndkf2*transpose(H)*inv((H*Pndkf2*transpose(H))+Rdkf2);
    xx2=xx2+Kx2*(z'-H*xx2);   
    xx_store_dkf2(:,i) = xx2(1:Nm);
    fd_est_dkf2(:,i) = fd_dkf2;
    % vv2=(xx2(1:Nn)-xxp2(1:Nn))/deltat;
    % aa2=(vv(1:Nn)-vvp2(1:Nn))/deltat;
    % FF=Mgm*aa(:)+Cgc*vv(:)+Kgk*xx(1:Nn);

    O_s2_st = obsv(Adkf2,H);
    
%% CHECKING THE OBSERVABILITY OF THE DAMAGES 
for z=1:Ndam
    Obsv_fd_s1(z) = Obsv_fd_s1(z) + norm(O_s1_fd(:, z))^2;
    
    Obsv_fd_s2(z) = Obsv_fd_s2(z) + norm(O_s2_fd(:, z))^2;
end

for k=1: 2*Nm

    Obsv_state_s1(k) = Obsv_state_s1(k) + norm(O_s1_st(:, k))^2;
    Obsv_state_s2(k) = Obsv_state_s2(k) + norm(O_s2_st(:, k))^2;
end
    
for k = 1: (2*Nm+Ndam)
    Obsv_ekf(k) = Obsv_ekf(k) + norm(O_ex(:,k))^2;

end
%     
%     

    
%    fuc(:,i)=inv(Hf)*(zf(:,i)-HBf);
    reba=1000;
    
   % fprintf(Fil1,'%f\n',fac');

      damfpredict(i,:)=fac;
   facc=fac;
%     
%     end
    P=(eye(dim)-K*H)*Pn;
    Pf=(eye(Ndam)-Kf*Hf)*Pnf;

     Pdkf2=(eye(dim)-Kx2*H)*Pndkf2;
    Pf2=(eye(Ndam)-Kdkf2*Hdkf2)*Pnf2;
   
    
xxp=xx;
vvp=vv;
aap=aa;

acp_2 = (xx2(Nm+1:2*Nm)-xxp2(Nm+1:2*Nm))/deltat;
xxp2 = xx2;


accp=acc;
   acp=ac;
   xp=x;
   vp=v;
   xrip=xri;
   vrip=vri;
   accpri=accri;

   
       
    
end
    
%% Vizualization 
% State plotting
figure(1);
sgtitle('Damage Estimation (sensor : z1,z2,z3,z4)');
for i=1:6
    subplot(2,3,i)
    plot(damf(i,:),'k-');
    hold on;
    plot(damfpredict(:,i),'b--'); hold on;
    plot(fd_est_dkf2(i,:),'m--');
    plot(fd_est_ex(i,:),'g--');
    xlabel('time');
    ylabel('Damage factor(fd)');
    hold off
end


% Damage factor plotting
figure(2);
for i=1:6
    subplot(2,3,i);
    plot(damf(i,:),'k-');
    hold on;
    plot(damfpredict(:,i),'b--'); hold on;
    plot(fd_est_dkf2(i,:),'m--');
    %plot(fd_est_ex(i,:),'g--');
    hold off;
end




% %% For the The observability matrix 
% 
% dt=0.00001;       % time step
% T_obs=0.1;    % observability window length (seconds)
% T = round(T_obs/dt);
% omega = sqrt(diag(Kgun)./diag(Mgm));
% F = Fgf*sin(w_1*(0:T-1)*dt);   % excite dominant mode
% 
% 
% % Initialize
% q  = zeros(Nm,N);
% qd = zeros(Nm,N);
% 
% 
% for j=1:Ndam
%     for i=1:Nm
%         for k=1:Nm
%             KKd(i,k)=Kd(j,i,k);
%         end
%     end
%     K_zone{j}=KKd;
% end
% omega = sqrt(diag(Kgun)./diag(Mgm));
% 
% Minv = inv(Mgm);
% F = Fgf*sin(w_1*(0:N-1)*deltat);
% for n = 2:N
%     qdd = Minv*(F(:,n-1) - Kgk*q(:,n-1));
%     qd(:,n) = qd(:,n-1) + deltat*qdd;
%     q(:,n)  = q(:,n-1)  + deltat*qd(:,n-1);
% end
% % For state obs
% for i = 1 : Ndam
%     G = eye(Nm)' * K_zone{i} * eye(Nm);      % Nm x Nm
%     H_at = Phi * Minv * G;         % nOberv x Nm
% 
%     for n = 1:N
%         D_i(:,n) = H_at * q(:,n);       % nOberv x 1
%     end
%     O(i) = mean(vecnorm(D_i,2,1));      % observability index
% end
% 
% for i = 1 : 2*Nm
%     % G = eye(Nm)' * K_zone{i} * eye(Nm);      % Nm x Nm
%     Hs = [Phi , X ; X , Phi] ;         % nOberv x Nm
% 
%     for n = 1:N
%         DS_i(:,n) = Hs * [q(:,n); qd(:,n) ];     % nOberv x 1
%     end
%     OS(i) = mean(vecnorm(DS_i,2,1));      % observability index
% end
% 
% OS = OS / max(OS);          % normalize to [0,1]
% epsiS = 1e-12;             % avoid divide-by-zero
% p10 = 1e-7;               % base variance
% P1_0 = diag(p10 ./ (OS + epsiS));
% q10 = 1e-8;               % base process noise
% Q1 = diag(q10 ./ (OS + epsiS));
% 
% P1_0 = min(P1_0, 1e2*eye(2*Nm));
% Q1   = min(Q1,   1e-1*eye(2*Nm));
% 
% 
% O = O / max(O);          % normalize to [0,1]
% epsi = 1e-6;             % avoid divide-by-zero
% 
% p0 = 1e-12;               % base variance
% P2_0 = diag(p0 ./ (O + epsi));
% q0 = 1e-15;               % base process noise
% Q2 = diag(q0 ./ (O + epsi));
% 
% P2_0 = min(P2_0, 1e2*eye(Ndam));
% Q2   = min(Q2,   1e-1*eye(Ndam));
% 
% figure(3);
% bar(OS);
% xlabel('No. of modes(disp & Vel) ');
% ylabel('Normalized observability');
% title('Zone-wise damage observability');
% grid on;
% 
% figure(4);
% bar(O);
% xlabel('Damage zone');
% ylabel('Normalized observability');
% title('Zone-wise damage observability');
% grid on;

% PLOTTING THE STATES  

figure;
  for i =1:Nm
        subplot(2,3,i);
        plot(linspace(0,deltat*N,N),xx_true_store(i,:),'k-'); hold on;
        plot(linspace(0,deltat*N,N),xx_store_dkf2(i,:),'m--'); 
        plot(linspace(0,deltat*N,N),xx_store_dkf(i,:),'b--'); 
        plot(linspace(0,deltat*N,N),xx_store_ekf(i,:),'g--');
        xlabel('time');
        ylabel('displacement');
        title('disp(damage) vs time');

  end


%% Plotting the observation matrix 
figure(6)
bar(Obsv_fd_s1);%/max(Obsv_fd_s1));
set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6'});
title('Observability(Fd) over Time for scheme1');

figure(7)
bar(Obsv_fd_s2);%/max(Obsv_fd_s2));
set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6'});
title('Observability(Fd) over Time for scheme2 ');
figure(8)
bar(Obsv_state_s1);%/max(Obsv_state_s1));
set(gca, 'XTickLabel', {'m 1', 'm 2', 'm 3', 'm 4', 'm 5', 'm 6', 'm 7','m 8','m9','m10','m11','m12'});
title(' Observability Strength (modes) over Time for scheme1 ');

figure(9)
bar(Obsv_state_s2);%/max(Obsv_state_s2));
set(gca, 'XTickLabel', {'m 1', 'm 2', 'm 3', 'm 4', 'm 5', 'm 6', 'm 7','m 8','m9','m10','m11','m12'});
title(' Observability Strength (modes) over Time for scheme2 ');


figure(10)
bar(Obsv_ekf(2*Nm+1:2*Nm+Ndam));%/max(Obsv_ekf(2*Nm+1:2*Nm+Ndam)));
set(gca, 'XTickLabel', { 'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6'});
title(' Observability Strength (dams) over Time for EKF ');


P2_s2 = 1./Obsv_fd_s2;
P2_s1 = 1./Obsv_fd_s1;

P1_s1 = 1./Obsv_state_s1;
P1_s2 = 1./Obsv_state_s2;


P_ekf = 1./Obsv_ekf ; 

%% printing 
for i =1:18
    P_ekf(i)
    %P2_s1(i)
end
