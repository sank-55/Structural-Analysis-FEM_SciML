%% Done modal domain with observability 


clear all
syms tau
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
w_1 = 20*pi; %1*137.2892;
w_6 =  897.6603; 

N = 5000;
deltat = 0.00005;
nObserv= 5 ;       % No of observable zones (** sensor is acting) 

syms x y

%Fil1=fopen('facdat1.m','w');


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
 % sensox=[a/6 a/2 5*a/6  a/6  5*a/6   a/2 ];
 % sensoy=[b/4 b/4 b/4  3*b/4  3*b/4  b/2 ];


 % Effective for 5 sensor location 
 sensox=[  a/6   5*a/6   1*a/6   5*a/6    a/2   ];
 sensoy=[   b/4   b/4    3*b/4    3*b/4    b/2 ];

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

%Phi=eye(nObserv);     
% Phi = [ 1  0  0  0  0  0;
%         0  1  0  0  0  0 ;
%         0  0  1  0  0  0 ;
%         0  0  0  1  0  0 ;
%         0  0  0  0  1  0 ;
%         0  0  0  0  0  1 ; 
%           ];
% X=zeros(nObserv,Nm);

%This is 4 sensor
% Phi = [ 1 0  0   0  0  0;
%         0  1  0  0  0  0;
%         0  0  1  0  0  0;
%         0  0  0   0  1  0   ];
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

%Pf=eye(Ndam)*10000;  % BEST AT 10000 
%Pf = 1e-12*blkdiag( 0.4793  , 0.3464 ,0.5248, 0.4994 ,0.3564, 0.5304);% for 6 sensor scheme 1;
% Pf2 = 1e-12*blkdiag(0.4809 , 0.3474 ,0.5266 , 0.5011,  0.3574 , 0.5322); % for 6 sensor scheme 2


% --- 5 sensor based 6 Ndam --------- 
Pf = 1e-14*blkdiag(0.0690,0.1179 ,0.0865,0.0895, 0.1190,  0.0968); % for 4 sensor 
Pf2 = 1e-14*blkdiag(0.0693,0.1181, 0.0869,0.0899,  0.1194 , 0.0972);



%Pf=load('Pf.txt');
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
Pekf = blkdiag(1e8*P_state_ekf,Q_ekf_fd);

Qekf = 1e6*blkdiag(3.3898e-15,3.8276e-17, 2.4722e-16, 1.3198e-17, 1.8844e-18,2.3570e-17, 3.4714e-08 ,3.8733e-10, 2.5383e-09,  1.3562e-10, 1.8996e-11, 2.4181e-10, 1.9583e-13, 1.3443e-13 ,  2.0059e-13, 1.9658e-13, 1.3482e-13, 2.0093e-13  );


%Qf=    ;
%Qf = 1e-11*blkdiag( 0.2438  ,0.1030 ,0.2463  ,0.2365   , 0.1000, 0.2328);  % working better in 1e-12
Qf = Pf; %1e-12*blkdiag( 0.4793  , 0.3464 ,0.5248, 0.4994 ,0.3564, 0.5304);
Qf2 =Pf2; %    % 1e1*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);

% erR=z-[x;v];
% erf=zf-acc;
% R=cov(erR')*10000000;
% Rf=cov(erf')*10000000;
R=1*eye(2*nObserv)*10^(-10);  % for state ( disp & vel) of Scheme - 1 for 6 sen ->1*eye(2*nObserv)*10^(-4) works better ( no noise)  
Rf=1*eye(nObserv)*10^(-10);
Rf2=1*eye(nObserv)*10^(-10);
%Rf = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);

Rdkf2=1*eye(2*nObserv)*10^(-10); % for state ( disp & vel) of Scheme - 2
%Rf2 = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);

Rekf = 1e1*eye(nObserv);
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
xx_store_dkf = zeros(2*Nm,N);
xx_store_dkf2 = zeros(2*Nm,N);
xx_store_ekf = zeros(2*Nm,N);

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
    
    xrri=Phi*xri';
    vrri=Phi*vri';
    z = [xrri; vrri]';
   %z(1:2*nObserv)=[xri(1:nObserv)'; vri(1:nObserv)'];
  % z(1:2*nObserv,i)=[ v(1:nObserv,i) ; ac(1:nObserv,i)];
   zf(1:nObserv)=Phi*ac;

   





%for i=2:N
    
    iii=i;
    t=(i-1)*deltat;
    
   
%     err=10000;
%     
%     
%     while err>0.001
%    
%       

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

% II=eye(Nn);
% 
% AA=Mgm+gamma*Cgc+beta*deltat^2*Kgk;
% AAinv=inv(AA);
% CC=Cgc+Kgk*deltat;
% BB=deltat*Cgc*(1-gamma)+deltat^2/2*(1-2*beta)*Kgk;
% 
% A11=zeros(3*Nm,3*Nm);
% 
% A11(1:Nn,:)=[II-beta*deltat^2*AAinv*Kgk deltat*II-beta*deltat^2*AAinv*CC (1-2*beta)/2*deltat^2*II-beta*deltat^2*AAinv*BB];
% 
% A11(Nn+1:2*Nn,:)=[-gamma*deltat*AAinv*Kgk II-gamma*deltat*AAinv*CC (1-gamma)*deltat*II-gamma*deltat*AAinv*BB];
% 
% A11(2*Nn+1:3*Nn,:)=AAinv*[-Kgk -CC -BB];
% A=A11;
% 
% B=[eye(Nm,Nm); gamma*deltat*eye(Nm,Nm); beta*deltat^2*eye(Nm,Nm)];

% A=exp(Ac*deltat);
% 
% B=double(int(exp(Ac*tau),tau,0,deltat))*Bc;

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
sgtitle('Damage Estimation (sensor : z1,z3,z4,z6,Mid)');
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




%% For the The observability matrix 

dt=0.00001;       % time step
T_obs=0.1;    % observability window length (seconds)
T = round(T_obs/dt);
omega = sqrt(diag(Kgun)./diag(Mgm));
F = Fgf*sin(w_1*(0:T-1)*dt);   % excite dominant mode


% Initialize
q  = zeros(Nm,N);
qd = zeros(Nm,N);


for j=1:Ndam
    for i=1:Nm
        for k=1:Nm
            KKd(i,k)=Kd(j,i,k);
        end
    end
    K_zone{j}=KKd;
end
omega = sqrt(diag(Kgun)./diag(Mgm));

Minv = inv(Mgm);
F = Fgf*sin(w_1*(0:N-1)*deltat);
for n = 2:N
    qdd = Minv*(F(:,n-1) - Kgk*q(:,n-1));
    qd(:,n) = qd(:,n-1) + deltat*qdd;
    q(:,n)  = q(:,n-1)  + deltat*qd(:,n-1);
end
% For state obs
for i = 1 : Ndam
    G = eye(Nm)' * K_zone{i} * eye(Nm);      % Nm x Nm
    H_at = Phi * Minv * G;         % nOberv x Nm

    for n = 1:N
        D_i(:,n) = H_at * q(:,n);       % nOberv x 1
    end
    O(i) = mean(vecnorm(D_i,2,1));      % observability index
end

for i = 1 : 2*Nm
    % G = eye(Nm)' * K_zone{i} * eye(Nm);      % Nm x Nm
    Hs = [Phi , X ; X , Phi] ;         % nOberv x Nm

    for n = 1:N
        DS_i(:,n) = Hs * [q(:,n); qd(:,n) ];     % nOberv x 1
    end
    OS(i) = mean(vecnorm(DS_i,2,1));      % observability index
end

OS = OS / max(OS);          % normalize to [0,1]
epsiS = 1e-12;             % avoid divide-by-zero
p10 = 1e-7;               % base variance
P1_0 = diag(p10 ./ (OS + epsiS));
q10 = 1e-8;               % base process noise
Q1 = diag(q10 ./ (OS + epsiS));

P1_0 = min(P1_0, 1e2*eye(2*Nm));
Q1   = min(Q1,   1e-1*eye(2*Nm));


O = O / max(O);          % normalize to [0,1]
epsi = 1e-6;             % avoid divide-by-zero

p0 = 1e-12;               % base variance
P2_0 = diag(p0 ./ (O + epsi));
q0 = 1e-15;               % base process noise
Q2 = diag(q0 ./ (O + epsi));

P2_0 = min(P2_0, 1e2*eye(Ndam));
Q2   = min(Q2,   1e-1*eye(Ndam));

figure(3);
bar(OS);
xlabel('No. of modes(disp & Vel) ');
ylabel('Normalized observability');
title('Zone-wise damage observability');
grid on;

figure(4);
bar(O);
xlabel('Damage zone');
ylabel('Normalized observability');
title('Zone-wise damage observability');
grid on;

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

%% printing 
for i =1:18
    P_ekf(i)
    %P2_s1(i)
end
