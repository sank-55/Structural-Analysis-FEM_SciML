% %% Initial 
% clear all; close all;
% syms tau
% p=1000;
% ro=2700; % ( in kg/m3)
% yo =70*10^9;%youngs modulus
% neu=0.3;
% h=0.001;%thickness
% a=0.6; %( m)
% b=0.3; % (m)
% divx=3;%input(' enter the no of divisions on the breadth ');
% % ex=length of element in x direction
% ex=a/(divx);
% divy=3;%input(' enter the no of divisions on the length ');
% % ey=length of element in y direction
% ey=b/(divy);
% %% natural frequencies at those given params
% w_1 = 45*pi ; %137.2892;
% w_6 =    20*pi;% 2*w_1; %897.6603; % 6th natural freq for this system 
% w_8 = 8*w_1;
% 
% N=500;
% deltat=0.00005;
% nObserv=5;
% 
% syms x y
% 
% %Fil1=fopen('facdat1.m','w');
% 
% 
% Nm=6;
% Nn=Nm;
% %phi=[sin(pi*x/a)*sin(pi*y/b) sin(3*pi*x/a)*sin(pi*y/b) sinHafpi*x/a)*sin(3*pi*y/b) sin(3*pi*x/a)*sin(3*pi*y/b) sin(4*pi*x/a)*sin(4*pi*y/b)];
% % For every Phi in serial manner 
% phi=[sin(pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(2*pi*y/b) sin(2*pi*x/a)*sin(pi*y/b) sin(2*pi*x/a)*sin(2*pi*y/b) sin(3*pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(3*pi*y/b)]; 
% % Onlytking the odd terms because evenmuliples integrated becomes zero
% %phi=[sin(pi*x/a)*sin(pi*y/b)  sin(pi*x/a)*sin(3*pi*y/b) sin(3*pi*x/a)*sin(pi*y/b)  sin(3*pi*x/a)*sin(3*pi*y/b) sin(pi*x/a)*sin(5*pi*y/b) sin(5*pi*x/a)*sin(1*pi*y/b) ]; 
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
% % Damage factor applying
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
% Kgun=Kgk;  % kGK Is the undamaged one 
% Haf=[1 0 0 0];
% 
% damf=zeros(Ndam,N);
% % 
%   for i=20:N
%       damf(1,i)=(i-20)^.7*0.3/(N-20)^.7;
%   end
% % % 
%   for i=120:N
%       damf(2,i)=(i-120)*0.6/(N-120);
%   end
% %  
%   for i=60:N
%      damf(3,i)=(i-60)*0.3/(N-60);
%  end
% % 
% for i=1:N
%     damf(4,i)=(i-1)*0.3/(N-1);
% end
% % 
% for i=300:N
%     damf(5,i)=(i-300)*0.2/(N-300);
% end
% % 
% 
% for i=3:N
%     damf(6,i)=(i-3)*0.2/(N-3);
% end
% 
% 
% dim=2*Nm;
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
%    % for serial phi including even terms 
%    %Phi(i,:)=[sin(pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(2*pi*y1/b) sin(2*pi*x1/a)*sin(pi*y1/b) sin(2*pi*x1/a)*sin(2*pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b)];
%     Phi(i,:)=[sin(pi*x1/a)*sin(pi*y1/b)  sin(pi*x1/a)*sin(3*pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b)  sin(3*pi*x1/a)*sin(3*pi*y1/b) sin(1*pi*x1/a)*sin(5*pi*y1/b) sin(5*pi*x1/a)*sin(1*pi*y1/b)];
% end
% 
% %Dekf constant Initializaion
% %Phi=eye(Nm);
% 
% % Phi = [ 1, 0, 0, 0, 0, 0 ;
% %         0, 1, 0, 0, 0, 0;
% %         0, 0, 0, 1, 0, 0;
% %         0, 0, 0, 0, 1, 0;
% %         0, 0, 0, 0, 0, 1 ];
% 
% X=zeros(nObserv,Nm);
% 
% % H=[Phi X X;X Phi X;X X Phi];
% 
% H=[Phi X;X Phi];
% 
% Haf=Phi;
% P=[ones(Nm,1)*10^(-2);ones(Nm,1)*10^(-1)];
% P=diag(P);
% Pdkf2=P;
% %Pf=eye(Ndam)*10000;  % BEST AT 10000 
% % Pf = 1e6*blkdiag(0.0100,0.0283,0.0714 ,0.0192,0.0230, 0.0209);
% % Pf2 = 1e6*blkdiag(0.0100,0.0283,0.0714 ,0.0192,0.0230, 0.0209);
% Pf = 1e-3*blkdiag(0.0100,0.0220,0.0100 ,0.0100,0.0220, 0.0100);
% Pf2 = Pf ;%1e6*blkdiag(0.0100,0.0283,0.0714 ,0.0192,0.0230, 0.0209);
% %Pf=load('Pf.txt');
% Q=eye(dim)/10^-8;  % 10^8
% 
% Qdkf2=eye(dim)/10^-5; 
% %Q(2*Nn+1:dim,2*Nn+1:dim)=eye(dim-(2*Nn));
% dimf=1;
% 
% Pekf = 1e-2*eye(2*Nm + Ndam);
% Qekf = 1e-5*eye(2*Nm + Ndam);
% 
% 
% %Qf=load('Qf.txt');
% % Qf = 1e-4*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);  % working better in 1e-2
% % Qf2 = 1e-4*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);
% Qf = 8e-5*blkdiag(0.0999,0.2197 , 0.0999 ,0.0999,0.2197,0.0999);
% Qf2 = Qf;  % 1e-4*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);
% % erR=z-[x;v];
% % erf=zf-acc;
% % R=cov(erR')*10000000;
% % Rf=cov(erf')*10000000;
% R=eye(2*nObserv)*10^(-1);
% Rf=eye(nObserv)*10^(11);
% Rdkf2=eye(2*nObserv)*10^(-1);
% Rf2=eye(nObserv)*10^(-1);
% 
% Rekf = 1e-1*eye(nObserv);
% 
% 
% %Dekf begins
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
% fd_ex = zeros(Ndam,1);
% fd_est_ex = zeros(Ndam,N);
% fd_dkf2 = zeros(Ndam,1);
% fd_est_dkf2 = zeros(Ndam,N);
% 
% xxp=zeros(2*Nm,1);
% xxp2=zeros(2*Nm,1);
% xxp_exfd=zeros(2*Nm+Ndam,1);  % for the extended kalman filter 
% vvp=zeros(Nm,1);
% aap=zeros(Nm,1);
% accp=zeros(Nm,1);
% acp=zeros(Nm,1);
% acp_2 = zeros(Nm,1);
% xp=zeros(Nm,1);
% vp=zeros(Nm,1);
% xrip=zeros(Nm,1);
% vrip=zeros(Nm,1);
% accpri=zeros(Nm,1);
% xx_store_dkf = zeros(2*Nm,N);
% xx_store_dkf2 = zeros(2*Nm,N);
% xx_store_ekf = zeros(2*Nm,N);
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
%     BB=Cgc*deltat*(1-gamma)+Kgk*(1-2*beta)*deltat^2/2;
% 
%     CC=Cgc+deltat*Kgk;
%     AAinv=inv(AA);
% 
%     t=(i-1)*deltat;
%     t;
%     Fg(:)=Fgf*sin(w_6*t);
%     acc(:)=AAinv*(Fg(:)-[BB]*accp(:)-[CC]*v(:)-Kgk*x(:));
%     accp=acc;
% 
%     for kk=1:Nm
%         ac(kk)=acc(kk)+ 0*(rand-.5)/50;
%         Accel(i,kk)=acc(kk);
%     end
%     accri=ac+ 0*(rand-.5)/70;
%     accri(1)=ac(1)+ 0*(rand-0.5)/5;
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
%     %z(1:2*nObserv)=[xri(1:nObserv)'; vri(1:nObserv)'];
%     z = [(Phi*xri'); (Phi*vri')]';
%    zf =( Phi * ac)';
% 
%     iii=i;
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
%     F=Fgf*sin(w_6*t);    
%     F1=F;
%     Mgminv=inv(Mgm);
% 
% 
%     Af=eye(Ndam);
%     Hf=Haf*Mgminv*Kdx;
%     Cgcx=Mgminv*(Kgun*disp+Cgc*vel);
%     HBf=Haf*(Mgminv*F-Cgcx);
% 
%     % for Scheme 2 
% 
%      Kgks2 = Kgun;
%      for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         Kgks2=Kgks2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
%      end
% 
%       Kdx2=zeros(Nm,Ndam);
% 
%     for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         xi(:,j)=tem*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2 );
%         Kdx2(:,j)=xi(:,j);
%     end
%     Hdkf2 = Phi*(inv(Mgm + deltat*0.5*Cgc + deltat^2*0.25*(Kgks2))*(Kdx2)) ;
%     a_dkf2 = Hdkf2*fd_dkf2 + Phi*(inv(Mgm + deltat*0.5*Cgc + deltat^2*0.25*(Kgks2))*(F1- Kgun*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2) - Cgc*(xxp2(Nm+1:2*Nm) + deltat^2*0.5*acp_2) ));
% 
% 
%     % for scheme 1 
%     Pnf=Af*Pf*transpose(Af)+Qf;
%     Kf=Pnf*transpose(Hf)*inv((Hf*Pnf*transpose(Hf))+Rf);
% 
%     %Af=eye(1);
%         reba=100;
%         fac=facc+Kf*(zf'-Hf*fac(:)-HBf);
% 
%      %for Scheme 2
%      Pnf2=Af*Pf2*transpose(Af)+Qf2;
%      Kdkf2=Pnf2*transpose(Hdkf2)*inv((Hdkf2*Pnf2*transpose(Hdkf2))+Rf2);
% 
%     %Af=eye(1);
%     reba=100;
%     fd_dkf2=fd_dkf2+Kdkf2*(zf'-a_dkf2);
% 
% %% 2nd Loop for DEKF 
% 
%     Kgk1 = Kgun;
%      for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         Kgk1=Kgk1-fac(j)*tem;  % here must be fac or facc damf will not be there
%      end
%      Ac=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgk1  -inv(Mgm)*Cgc];
%      Bc=[zeros(Nm,Nm); inv(Mgm)];
% 
%  % FOr Scheme 2 
% 
%     Kgdkf2 = Kgun;
%      for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         Kgdkf2=Kgdkf2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
%      end
% 
% 
%     Adkf2=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgdkf2  -inv(Mgm)*Cgc];
%      Bdkf2=[zeros(Nm,Nm); inv(Mgm)];
% 
% 
% 
%     % Extended kalman filter matrix 
%      S_ekf=zeros(Nm,Ndam);
%     for j=1:Ndam
%        for k=1:Nm
%            for m=1:Nm
%                tem(k,m)=Kd(j,k,m);
%            end
%        end
%        S_ekf(:,j) = tem*xxp_exfd(1:Nm);
%     end
% 
%     Aex = [zeros(Nm,Nm) , eye(Nm,Nm) , zeros(Nm,Ndam);
%             -inv(Mgm)*Kgk1 , -inv(Mgm)*Cgc , inv(Mgm)*S_ekf ;
%             zeros(Ndam,Nm), zeros(Ndam,Nm), eye(Ndam,Ndam)];
% 
%     Bex = [zeros(Nm,Nm); inv(Mgm); zeros(Ndam,Nm)];
% 
% 
%         % B=[eye(Nm,Nm); gamma*deltat*eye(Nm,Nm); beta*deltat^2*eye(Nm,Nm)];
% 
%     % A=exp(Ac*deltat);
%     % 
%     % B=double(int(exp(Ac*tau),tau,0,deltat))*Bc;
% 
%     % Discretization 
%     % for scheme 1
%     A=deltat*Ac+eye(2*Nm);
%     B=deltat*Bc;
% 
%     % How this is defined
%     A=zeros(2*Nm,2*Nm);
%     BB=zeros(2*Nm,2*Nm);
%     B=zeros(2*Nm,Nm);
%     for jj=0:10
%         A=A+(Ac*deltat)^jj/factorial(jj);
%         BB=BB+(Ac*deltat)^jj/factorial(jj+1);
%     end
% 
%     B=BB*Bc*deltat;
% 
%     % for scheme 2 
%     Ax2=deltat*Adkf2+eye(2*Nm);
%     Bx2=deltat*Bdkf2;
% 
%     % How this is defined
%     Ax2=zeros(2*Nm,2*Nm);
%     Bx2_=zeros(2*Nm,2*Nm);
%     Bx2=zeros(2*Nm,Nm);
%     for jj=0:10
%         Ax2=Ax2+(Adkf2*deltat)^jj/factorial(jj);
%         Bx2_=Bx2_+(Adkf2*deltat)^jj/factorial(jj+1);
%     end
% 
%     Bx2=Bx2_*Bdkf2*deltat;
% 
% 
%     % For extended kalman filter 
%     Aex =  deltat*Aex+eye(2*Nm + Ndam);
%     Bex=deltat*Bex;
% 
%     Ax=zeros(2*Nm+Ndam,2*Nm+Ndam);
%     BBx=zeros(2*Nm+Ndam,2*Nm+Ndam);
%     Bx=zeros(2*Nm+Ndam,Nm);
%     for jj=0:10
%         Ax=Ax+(Aex*deltat)^jj/factorial(jj);
%         BBx=BBx+(Aex*deltat)^jj/factorial(jj+1);
%     end
% 
%     Bx=BBx*Bex*deltat;
% 
%     % for scheme 2 
%     Ax2=deltat*Adkf2+eye(2*Nm);
%     Bx2=deltat*Bdkf2;
% 
%     % How this is defined
%     Ax2=zeros(2*Nm,2*Nm);
%     Bx2_=zeros(2*Nm,2*Nm);
%     Bx2=zeros(2*Nm,Nm);
%     for jj=0:10
%         Ax2=Ax2+(Adkf2*deltat)^jj/factorial(jj);
%         Bx2_=Bx2_+(Adkf2*deltat)^jj/factorial(jj+1);
%     end
% 
%     Bx2=Bx2_*Bdkf2*deltat;
% 
% 
% % For extended kalman filter 
%     Aex =  deltat*Aex+eye(2*Nm + Ndam);
%     Bex=deltat*Bex;
% 
%     Ax=zeros(2*Nm+Ndam,2*Nm+Ndam);
%     BBx=zeros(2*Nm+Ndam,2*Nm+Ndam);
%     Bx=zeros(2*Nm+Ndam,Nm);
%     for jj=0:10
%         Ax=Ax+(Aex*deltat)^jj/factorial(jj);
%         BBx=BBx+(Aex*deltat)^jj/factorial(jj+1);
%     end
% 
%     Bx=BBx*Bex*deltat;
%      xx_exfd = Ax*xxp_exfd + Bx*F;
%     Pekfn = Ax*Pekf*transpose(Ax)+ Qekf; 
%      Kgk2 = Kgun;
% 
%      for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         Kgk2=Kgk2-xx_exfd(2*Nm+j)*tem;  % here must be fac or facc damf will not be there
%     end
%       S_ekf=zeros(Nm,Ndam);
%     for j=1:Ndam
%        for k=1:Nm
%            for m=1:Nm
%                tem(k,m)=Kd(j,k,m);
%            end
%        end
%        S_ekf(:,j) = tem*xx_exfd(1:Nm);
%     end
% 
% 
% 
%     Hekf = [-Phi*inv(Mgm)*Kgk2 , -Phi*inv(Mgm)*Cgc , Phi*inv(Mgm)*S_ekf ];
%     z_diff_ekf = zf' - Phi*( Mgm\(F1 - Kgk2*xx_exfd(1:Nm) - Cgc*xx_exfd(Nm+1 :2*Nm) )) ;
%     Kekf = Pekfn*transpose(Hekf)*inv((Hekf*Pekfn*transpose(Hekf)) + Rekf);
%     xxp_exfd = xx_exfd + Kekf*z_diff_ekf;
% 
%     % Storing the Values of the updated ekf 
%     xx_store_ekf = xxp_exfd(1:2*Nm); %(  first Nm 'u' then Nm 'vel')
%     fd_est_ex(:,i) = xxp_exfd(2*Nm +1: 2*Nm + Ndam);
% 
% 
%    % 0*B*Fg(:,i-1)
%    % For Scheme 1 
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
%     % For Scheme 2 
%     xx2=Ax2*xxp2+Bx2*F;
%     Pndkf2=Ax2*Pdkf2*transpose(Ax2)+Qdkf2;    
%     Kx2=Pndkf2*transpose(H)*inv((H*Pndkf2*transpose(H))+Rdkf2);
%     xx2=xx2+Kx2*(z'-H*xx2);   
%     xx_store_dkf2(:,i) = xx2;
%     fd_est_dkf2(:,i) = fd_dkf2;
%     % vv2=(xx2(1:Nn)-xxp2(1:Nn))/deltat;
% 
% 
%  reba=1000;
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
%      Pdkf2=(eye(dim)-Kx2*H)*Pndkf2;
%     Pf2=(eye(Ndam)-Kdkf2*Hdkf2)*Pnf2;
% 
% 
%     xxp=xx;
%     vvp=vv;
%     aap=aa;
% 
%     acp_2 = (xx2(Nm+1:2*Nm)-xxp2(Nm+1:2*Nm))/deltat;
%     xxp2 = xx2;
% 
% 
%     accp=acc;
%     acp=ac;
%     xp=x;
%     vp=v;
%     xrip=xri;
%     vrip=vri;
%     accpri=accri;
% 
% end
% 
% %% 
% for i=1:6
%     subplot(2,3,i)
%     plot(damf(i,:),'k-');
%     hold on;
%     plot(damfpredict(:,i),'b--'); hold on;
%     plot(fd_est_dkf2(i,:),'m--');
%     plot(fd_est_ex(i,:),'g--');
%     hold off
% end
% 




% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %    %% 12 Modes, 12 Damage Factors - Dual KF & EKF with Observability
% %
% %
% %%%%% _____________________________________________________________________________________________________________________________________________________________________________________________
% 
clear all; close all;
syms  x y

% 1. Physical Parameters
p = 1000;
ro = 2700;           % Density (kg/m3)
yo = 70*10^9;        % Young's modulus
neu = 0.3;           % Poisson's ratio
h = 0.001;           % Thickness (m)
a = 0.6;             % Length (m)
b = 0.4;             % Breadth (m)
D = (yo * h^3) / (12 * (1 - neu^2));

%% 2. Simulation Setup
w_1 = 100*pi;   %1 * 137.2892;
N = 5000;
deltat = 0.0001;
nObserv = 12;         % No of observable zones (sensors)

Nm = 6;             % 6 Modes
Nn = Nm;
Ndam = 12;           % 12 Damage zones

% 12 Mode Shapes Definition
phi = [sin(1*pi*x/a)*sin(1*pi*y/b), sin(1*pi*x/a)*sin(2*pi*y/b), ...
       sin(2*pi*x/a)*sin(1*pi*y/b), sin(2*pi*x/a)*sin(2*pi*y/b), ...
       sin(1*pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(1*pi*y/b)];
       % sin(2*pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(2*pi*y/b), ...
       % sin(3*pi*x/a)*sin(3*pi*y/b), sin(1*pi*x/a)*sin(4*pi*y/b), ...
       % sin(4*pi*x/a)*sin(1*pi*y/b), sin(2*pi*x/a)*sin(4*pi*y/b)];

dxxphi = diff(phi, x, 2);
dxyphi = diff(diff(phi, x), y);
dyyphi = diff(phi, y, 2);

% 12 Damage Zone Coordinates (4x3 Grid over the plate)
damdox = [0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a; ...
          0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a; ...
          0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a];

damdoy = [0 b/3; 0 b/3; 0 b/3; 0 b/3; ...
          b/3 2*b/3; b/3 2*b/3; b/3 2*b/3; b/3 2*b/3; ...
          2*b/3 b; 2*b/3 b; 2*b/3 b; 2*b/3 b];

% 4 Sensor Coordinates
% %             z1      z2          z3      z4        z5    z6       z7       z8       z9      z10     z11      z12     BM      MM        TM
% sensox = [  a/8      3*a/8     5*a/8    7*a/8     a/8    3*a/8    5*a/8    7*a/8    a/8             5*a/8                     a/2       a/2        ];%    a/2                 ];
% sensoy = [  b/6       b/6        b/6      b/6     b/2      b/2     b/2       b/2    5*b/6            5*b/6                     b/2       5*b/6      ];%   b/6            ];
% 

%% sensor from the highest observibility 
sensox = [ 0.1220   0.4780   0.1220  0.4780  0.4881  0.1119  0.2949  0.3051  0.3966  0.2034   0.3254  0.2136 ];
sensoy = [ 0.0718  0.0718  0.3282   0.3282  0.2359  0.2359  0.3282   0.0718  0.1538   0.1538   0.2359  0.2564 ];

%% 3. Matrix Pre-allocation & Symbolic Integration
Mgm = zeros(Nm, Nm);
Kgk = zeros(Nm, Nm);
Fgf = zeros(Nm, 1);
Kd = zeros(Ndam, Nm, Nm);

for i = 1:Nm
    for j = 1:Nm
        Mgm(i,j) = double(ro * h * int(int(phi(i)*phi(j), x, 0, a), y, 0, b));
        Kgk(i,j) = double(D * int(int(dxxphi(i)*dxxphi(j) + 2*dxyphi(i)*dxyphi(j) + dyyphi(i)*dyyphi(j), x, 0, a), y, 0, b));
    end
    Fgf(i,1) = double(int(int(p*phi(i), x, 0, a), y, 0, b));
end

Cgc = 0.0003 * Mgm + 0.0003 * Kgk;     

for i = 1:Ndam
    for j = 1:Nm
        for k = 1:Nm
            Kd(i,j,k) = double(D * int(int(dxxphi(j)*dxxphi(k) + 2*dxyphi(j)*dxyphi(k) + dyyphi(j)*dyyphi(k), x, damdox(i,1), damdox(i,2)), y, damdoy(i,1), damdoy(i,2)));
        end
    end
end

%% 4. Initial Conditions & Damage Scenarios
gamma = 1/2;
beta = 1/2;
x_st = zeros(Nm,1);
v_st = zeros(Nm,1);
acc = zeros(Nm,1);
ac = zeros(Nm,1);
Fg = zeros(Nm,1);
Kgun = Kgk;  

% Damage Injection (Spanning up to 12 zones)
damf = zeros(Ndam, N);
for i = 200:N,  damf(1,i) = (i-200)^.7 * 0.3 / (N-200)^.7; end
for i = 1200:N, damf(2,i) = (i-1200) * 0.6 / (N-1200); end
for i = 600:N,  damf(3,i) = (i-600) * 0.3 / (N-600); end
for i = 30:N,   damf(5,i) = (i-30) * 0.2 / (N-30); end
for i = 800:N,  damf(8,i) = (i-800) * 0.4 / (N-800); end
for i = 1500:N, damf(11,i) = (i-1500) * 0.25 / (N-1500); end

%% 5. Sensor Matrix Configuration (12 Modes mapped)
dim = 2 * Nm;
Phi = zeros(nObserv, Nm);
for i = 1:nObserv
    x1 = sensox(i);
    y1 = sensoy(i);
    Phi(i,:) = [sin(1*pi*x1/a)*sin(1*pi*y1/b), sin(1*pi*x1/a)*sin(2*pi*y1/b), ...
                sin(2*pi*x1/a)*sin(1*pi*y1/b), sin(2*pi*x1/a)*sin(2*pi*y1/b), ...
                sin(1*pi*x1/a)*sin(3*pi*y1/b), sin(3*pi*x1/a)*sin(1*pi*y1/b), ...
               ];
end

X = zeros(nObserv, Nm);
%Phi = eye(nObserv);
H = [Phi X ; X Phi];
Haf = Phi;

%% 6. Filter Tuning Matrices (Generalized for 12 Modes/Zones)
% Note: Multipliers retained from your code, structured for 12 dimensions
P = 1e-4 * eye(2*Nm);   
Pdkf2 = 1e-4 * eye(2*Nm); 

% Pf = 1e-8 * blkdiag(0.1909, 0.0999, 0.1044, 0.2145,0.1462 ,0.0498, 0.0506,0.1525 ,0.1894,0.1165, 0.1223,0.2172); 
% Pf2 = 1e-8 * blkdiag(0.1918, 0.1009, 0.1054,0.2156, 0.1472 , 0.0504, 0.0512,0.1536, 0.1902,0.1176,0.1235 ,0.2182);

Pf = 1e-8*blkdiag(0.1288,0.0468,0.0494,0.1424,0.1454, 0.0150,0.0151,0.1547,0.1187,0.0621,0.0638,0.1456);
Pf2 = 1e-8*blkdiag(0.1301,0.04768,0.0502,0.1438,0.1484, 0.0154,0.0155,0.1577,0.1197,0.0631,0.0649,0.1470);

Q = eye(dim) / 1e-12;  
Qdkf2 = eye(dim) / 1e-12; 

%Pekf = 1e19 * eye(2*Nm + Ndam); 
Pekf = 1e-1*blkdiag(5.3489e-15,6.7229e-17,4.0723e-16, 2.5243e-17,4.9120e-18,4.5925e-17,6.0385e-09,7.5071e-11,4.5675e-10 ,2.7876e-11,5.4014e-12, 5.0964e-11,2.4343e-08,1.4394e-08,1.4358e-08,2.4396e-08,1.5747e-08,6.0527e-09,6.1104e-09,1.6249e-08,1.9990e-08,1.3034e-08 ,1.4240e-08,2.4590e-08 );

Qekf = 1e-20 * eye(2*Nm + Ndam);
Qf = Pf; 
Qf2 = Pf2; 

R = 10^1.5 * eye(2*nObserv);  
Rf = 10^0 * eye(nObserv);
Rdkf2 = 10^1.5 * eye(2*nObserv); 
Rf2 = 10^0 * eye(nObserv);
Rekf = 10^-0.7 * eye(nObserv);

%% 7. Variable Initialization for Time Loop
xx = zeros(dim, 1);
Fo(1) = 0;
FF = zeros(Nm, 1);
aa = zeros(Nm, 1);
vv = zeros(Nm, 1);
fac = zeros(Ndam, 1);
facc = fac;
fd_dkf2 = zeros(Ndam, 1);
fd_est_ex = zeros(Ndam, N);
fd_est_dkf2 = zeros(Ndam, N);
xxp = zeros(2*Nm, 1);
xxp2 = zeros(2*Nm, 1);
xxp_exfd = zeros(2*Nm+Ndam, 1);  
vvp = zeros(Nm, 1);
aap = zeros(Nm, 1);
accp = zeros(Nm, 1);
acp = zeros(Nm, 1);
acp_2 = zeros(Nm, 1);
xp = zeros(Nm, 1);
vp = zeros(Nm, 1);
xrip = zeros(Nm, 1);
vrip = zeros(Nm, 1);
accpri = zeros(Nm, 1);
xx_store_true = zeros(Nm,N);
xx_store_dkf1 = zeros(Nm,N);
xx_store_dkf2 = zeros(Nm, N);
xx_store_ekf = zeros(Nm, N);

Obsv_state_s1 = zeros(2*Nm, 1);
Obsv_fd_s2 = zeros(Ndam, 1);
Obsv_fd_s1 = zeros(Ndam, 1);
Obsv_state_s2 = zeros(2*Nm, 1);
Obsv_ekf = zeros(2*Nm + Ndam, 1);

tem = zeros(Nm, Nm);
xi = zeros(Nm, Ndam);
damfpredict = zeros(N, Ndam);

%% 8. Main Time Integration Loop
for i = 2:N
    % True System Update
    Kgk = Kgun;
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        Kgk = Kgk - damf(j,i) * tem;
    end

    AA = Mgm + Cgc*deltat*gamma + Kgk*beta*deltat^2;
    BB = Cgc*deltat*(1-gamma) + Kgk*(1-2*beta)*deltat^2/2;
    CC = Cgc + deltat*Kgk;

    t = (i-1)*deltat;
    Fg(:) = Fgf * sin(w_1*t);
    % Using \ instead of inv() for numerical stability
    acc(:) = AA \ (Fg(:) - BB*accp(:) - CC*v_st(:) - Kgk*x_st(:));
    accp = acc;

    for kk = 1:Nm
        ac(kk) = acc(kk) + 1e-4*(randn-.5)/50;
    end
    accri = ac + 1e-4*(randn-.5)/70;
    accri(1) = ac(1) + 1e-4*(rand-0.5)/5;

    v_st(:) = vp + deltat*(1-gamma)*accpri + deltat*gamma*accri;
    x_st(:) = xp + deltat*vp + (1-2*beta)/2*deltat^2*accpri + beta*deltat^2*accri;
    vp = v_st;
    xp = x_st;

    for j = 1:Nm
        vri(j) = vrip(j) + (acp(j)+ac(j))/2 * deltat;
        xri(j) = xrip(j) + (vrip(j)+vri(j))/2 * deltat;
    end

    xx_store_true(:,i)=xri;

    z = [(Phi*(xri+xri*0.02*rand)'); (Phi*(vri+vri*0.02*rand)')]';
    zf =( Phi * (ac + ac*0.03*rand))';

    %% Filter Prediction & Update
    disp = xxp(1:Nm) + xxp(Nm+1:2*Nm)*deltat;
    vel = xxp(Nm+1:2*Nm) + acp(1:Nm)*deltat;
    fac = facc;
    Kdx = zeros(Nm, Ndam);

    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        xi(:,j) = tem * disp;
        Kdx(:,j) = xi(:,j);
    end

    F = Fgf * sin(w_1*t);    
    F1 = F;

    Af = eye(Ndam);
    Hf = Haf * (Mgm \ Kdx);
    Cgcx = Mgm \ (Kgun*disp + Cgc*vel);
    HBf = Haf * ((Mgm \ F) - Cgcx);
    O_s1_fd = obsv(Af, Hf);

    % Scheme 2 Stiffness
    Kgks2 = Kgun;
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        Kgks2 = Kgks2 - fd_dkf2(j) * tem; 
    end

    Kdx2 = zeros(Nm, Ndam);
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        xi(:,j) = tem * (xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2);
        Kdx2(:,j) = xi(:,j);
    end

    InvMatS2 = Mgm + deltat*0.5*Cgc + deltat^2*0.25*Kgks2;
    Hdkf2 = Phi * (InvMatS2 \ Kdx2);
    a_dkf2 = Hdkf2 * fd_dkf2 + Phi * (InvMatS2 \ (F1 - Kgun*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2) - Cgc*(xxp2(Nm+1:2*Nm) + deltat^2*0.5*acp_2)));
    O_s2_fd = obsv(Af, Hdkf2);

    % Update Damage Factors (Scheme 1)
    Pnf = Af * Pf * Af' + Qf;
    Kf = Pnf * Hf' / ((Hf * Pnf * Hf') + Rf);
    fac = facc + Kf * (zf' - Hf*fac(:) - HBf);

    % Update Damage Factors (Scheme 2)
    Pnf2 = Af * Pf2 * Af' + Qf2;
    Kdkf2 = Pnf2 * Hdkf2' / ((Hdkf2 * Pnf2 * Hdkf2') + Rf2);
    fd_dkf2 = fd_dkf2 + Kdkf2 * (zf' - a_dkf2);

    %% 2nd Loop for DEKF
    Kgk1 = Kgun;
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        Kgk1 = Kgk1 - fac(j)*tem;  
    end

    Ac = [zeros(Nm,Nm) eye(Nm,Nm); -(Mgm \ Kgk1) -(Mgm \ Cgc)];
    Bc = [zeros(Nm,Nm); Mgm \ eye(Nm)];

    % Scheme 2 Matrix
    Kgdkf2 = Kgun;
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        Kgdkf2 = Kgdkf2 - fd_dkf2(j)*tem;  
    end

    Adkf2 = [zeros(Nm,Nm) eye(Nm,Nm); -(Mgm \ Kgdkf2) -(Mgm \ Cgc)];
    Bdkf2 = [zeros(Nm,Nm); Mgm \ eye(Nm)];

    % Extended Kalman Filter Matrix
    S_ekf = zeros(Nm, Ndam);
    for j = 1:Ndam
       for k = 1:Nm
           for m = 1:Nm
               tem(k,m) = Kd(j,k,m);
           end
       end
       S_ekf(:,j) = tem * xxp_exfd(1:Nm);
    end

    Aex = [zeros(Nm,Nm) , eye(Nm,Nm) , zeros(Nm,Ndam);
           -(Mgm \ Kgk1) , -(Mgm \ Cgc) , (Mgm \ S_ekf);
           zeros(Ndam,Nm), zeros(Ndam,Nm), eye(Ndam,Ndam)];
    Bex = [zeros(Nm,Nm); Mgm \ eye(Nm); zeros(Ndam,Nm)];

    % Discretization Scheme 1
    A = zeros(2*Nm, 2*Nm);
    BB = zeros(2*Nm, 2*Nm);
    for jj = 0:10
        A = A + (Ac*deltat)^jj / factorial(jj);
        BB = BB + (Ac*deltat)^jj / factorial(jj+1);
    end
    B = BB * Bc * deltat;

    % Discretization Scheme 2
    Ax2 = zeros(2*Nm, 2*Nm);
    Bx2_ = zeros(2*Nm, 2*Nm);
    for jj = 0:10
        Ax2 = Ax2 + (Adkf2*deltat)^jj / factorial(jj);
        Bx2_ = Bx2_ + (Adkf2*deltat)^jj / factorial(jj+1);
    end
    Bx2 = Bx2_ * Bdkf2 * deltat;

    % Discretization EKF
    Ax = zeros(2*Nm+Ndam, 2*Nm+Ndam);
    BBx = zeros(2*Nm+Ndam, 2*Nm+Ndam);
    for jj = 0:10
        Ax = Ax + (Aex*deltat)^jj / factorial(jj);
        BBx = BBx + (Aex*deltat)^jj / factorial(jj+1);
    end
    Bx = BBx * Bex * deltat;

    % EKF Update
    xx_exfd = Ax * xxp_exfd + Bx * F;
    Pekfn = Ax * Pekf * Ax' + Qekf; 

    Kgk2 = Kgun;
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        Kgk2 = Kgk2 - xx_exfd(2*Nm+j) * tem;  
    end

    for j = 1:Ndam
       for k = 1:Nm
           for m = 1:Nm
               tem(k,m) = Kd(j,k,m);
           end
       end
       S_ekf(:,j) = tem * xx_exfd(1:Nm);
    end

    Hekf = [-Phi*(Mgm \ Kgk2) , -Phi*(Mgm \ Cgc) , Phi*(Mgm \ S_ekf) ];
    z_diff_ekf = zf' - Phi*(Mgm \ (F1 - Kgk2*xx_exfd(1:Nm) - Cgc*xx_exfd(Nm+1 :2*Nm)));
    Kekf = Pekfn * Hekf' / ((Hekf * Pekfn * Hekf') + Rekf);
    xxp_exfd = xx_exfd + Kekf * z_diff_ekf;

    xx_store_ekf(:,i) = xxp_exfd(1:Nm); 
    fd_est_ex(:,i) = xxp_exfd(2*Nm +1: 2*Nm + Ndam);
    O_ex = obsv(Ax, Hekf);

    % Scheme 1 State Update
    xx = A*xxp + B*F;
    Pn = A*P*A' + Q;    
    K = Pn * H' / ((H*Pn*H') + R);
    xx = xx + K*(z' - H*xx);    
    xx_store_dkf1(:,i) = xx(1:Nm);
    vv = (xx(1:Nn) - xxp(1:Nn)) / deltat;
    aa = (vv(1:Nn) - vvp(1:Nn)) / deltat;
    O_s1_st = obsv(Ac, H);

    % Scheme 2 State Update
    xx2 = Ax2*xxp2 + Bx2*F;
    Pndkf2 = Ax2*Pdkf2*Ax2' + Qdkf2;    
    Kx2 = Pndkf2 * H' / ((H*Pndkf2*H') + Rdkf2);
    xx2 = xx2 + Kx2*(z' - H*xx2);   
    xx_store_dkf2(:,i) = xx2(1:Nm);
    fd_est_dkf2(:,i) = fd_dkf2;
    O_s2_st = obsv(Adkf2, H);

    %% Accumulating Observability Metrics
    for z_idx = 1:Ndam
        Obsv_fd_s1(z_idx) = Obsv_fd_s1(z_idx) + norm(O_s1_fd(:, z_idx))^2;
        Obsv_fd_s2(z_idx) = Obsv_fd_s2(z_idx) + norm(O_s2_fd(:, z_idx))^2;
    end
    for k = 1:2*Nm
        Obsv_state_s1(k) = Obsv_state_s1(k) + norm(O_s1_st(:, k))^2;
        Obsv_state_s2(k) = Obsv_state_s2(k) + norm(O_s2_st(:, k))^2;
    end
    for k = 1:(2*Nm+Ndam)
        Obsv_ekf(k) = Obsv_ekf(k) + norm(O_ex(:,k))^2;
    end

    damfpredict(i,:) = fac;
    facc = fac;

    % Covariance Updates
    P = (eye(dim) - K*H) * Pn;
    Pf = (eye(Ndam) - Kf*Hf) * Pnf;
    Pdkf2 = (eye(dim) - Kx2*H) * Pndkf2;
    Pf2 = (eye(Ndam) - Kdkf2*Hdkf2) * Pnf2;

    xxp = xx;
    vvp = vv;
    aap = aa;
    acp_2 = (xx2(Nm+1:2*Nm) - xxp2(Nm+1:2*Nm)) / deltat;
    xxp2 = xx2;
    accp = acc;
    acp = ac;
    xrip = xri;
    vrip = vri;
    accpri = accri;
end

%% 9. Visualization & Plotting (Updated for 12 Damages)
figure(1);
sgtitle('Damage Factors: True vs Estimated (Zones 1-12)');
for i = 1:Ndam
    subplot(4,3,i)

    plot(damf(i,:), 'k-'); hold on;
    plot(damfpredict(:,i), 'b--');
    plot(fd_est_dkf2(i,:), 'm--');
    plot(fd_est_ex(i,:), 'g--');
    title(sprintf('Zone %d', i));
    xlabel('time');
    ylabel('Damage factor');
    hold off
end

%% 10. Secondary Observability Analysis
dt = 0.00001;       
T_obs = 0.1;    
T = round(T_obs/dt);
q = zeros(Nm,N);
qd = zeros(Nm,N);

for j = 1:Ndam
    for i = 1:Nm
        for k = 1:Nm
            KKd(i,k) = Kd(j,i,k);
        end
    end
    K_zone{j} = KKd;
end

F = Fgf * sin(w_1*(0:N-1)*deltat);

for n = 2:N
    qdd = Mgm \ (F(:,n-1) - Kgk*q(:,n-1));
    qd(:,n) = qd(:,n-1) + deltat*qdd;
    q(:,n) = q(:,n-1) + deltat*qd(:,n-1);
end

for i = 1:Ndam
    G = eye(Nm)' * K_zone{i} * eye(Nm);      
    H_at = Phi * (Mgm \ G);         
    for n = 1:N
        D_i(:,n) = H_at * q(:,n);       
    end
    O(i) = mean(vecnorm(D_i,2,1));      
end

for i = 1:2*Nm
    Hs = [Phi , X ; X , Phi];         
    for n = 1:N
        DS_i(:,n) = Hs * [q(:,n); qd(:,n)];     
    end
    OS(i) = mean(vecnorm(DS_i,2,1));      
end

OS = OS / max(OS);          
O = O / max(O);          

%% 11. Final Plotting (Expanded for 12/24 ranges)
figure(2);
bar(OS);
xlabel('No. of modes (disp & Vel: 1-24)');
ylabel('Normalized observability');
title('State Observability');
grid on;

figure(3);
bar(O);
xlabel('Damage zone (1-12)');
ylabel('Normalized observability');
title('Zone-wise damage observability');
grid on;

% Label Generics for 12/24 sizes
zone_labels = arrayfun(@(x) sprintf('Z%d', x), 1:12, 'UniformOutput', false);
mode_labels = arrayfun(@(x) sprintf('m%d', x), 1:24, 'UniformOutput', false);

figure(4)
bar(Obsv_fd_s1);
set(gca, 'XTickLabel', zone_labels);
title('Observability (Fd) over Time for scheme 1');

figure(5)
bar(Obsv_fd_s2);
set(gca, 'XTickLabel', zone_labels);
title('Observability (Fd) over Time for scheme 2');

figure(6)
bar(Obsv_state_s1);
set(gca, 'XTickLabel', mode_labels);
title('Observability Strength (modes) over Time for scheme 1');

figure(7)
bar(Obsv_state_s2);
set(gca, 'XTickLabel', mode_labels);
title('Observability Strength (modes) over Time for scheme 2');

figure(8)
bar(Obsv_ekf(2*Nm+1:2*Nm+Ndam));
set(gca, 'XTickLabel', zone_labels);
title('Observability Strength (dams) over Time for EKF');

P1_s1 = 1./Obsv_state_s1;
P1_s2 = 1./Obsv_state_s2;

P2_s1 = 1./Obsv_fd_s1;
P2_s2 = 1./Obsv_fd_s2;
%%
figure(10);
  for i =1:Nm
        subplot(2,3,i);
        plot(linspace(0,deltat*N,N),xx_store_true(i,:),'k-'); hold on;
        plot(linspace(0,deltat*N,N),xx_store_dkf2(i,:),'m--'); hold on; 
        plot(linspace(0,deltat*N,N),xx_store_dkf1(i,:),'b--'); hold on;
        plot(linspace(0,deltat*N,N),xx_store_ekf(i,:),'g--');
        xlabel('time');
        ylabel('displacement');
        title('displacement vs time');

  end


%% 
% P_ekf = 1./Obsv_ekf;
% for i=1:24
%     P_ekf(i)
% end

%% cHECKING THE REDUNDANCY 
% [R, pivots] = rref(Phi);
% 
% % 'pivots' will list the columns that are independent.
% % Any number from 1 to 12 NOT in 'pivots' is a redundant column.
% redundant_cols = setdiff(1:12, R);
% disp('Redundant Columns are:');
% disp(redundant_cols);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   1000 load and 100 hz freq 12 s 12 zones 6 mode shapes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %    %% 12 Modes, 12 Damage Factors - Dual KF & EKF with Observability
% %
% %
% %%%%% _____________________________________________________________________________________________________________________________________________________________________________________________
% 
clear all; close all;
syms  x y

% 1. Physical Parameters
p = 1000;
ro = 2700;           % Density (kg/m3)
yo = 70*10^9;        % Young's modulus
neu = 0.3;           % Poisson's ratio
h = 0.001;           % Thickness (m)
a = 0.6;             % Length (m)
b = 0.4;             % Breadth (m)
D = (yo * h^3) / (12 * (1 - neu^2));

%% 2. Simulation Setup
w_1 = 100*pi;   %1 * 137.2892;
N = 5000;
deltat = 0.0001;
nObserv = 12;     %(8 to 12 easily predicted by it  )     % No of observable zones (sensors)

Nm = 6;             % 12 Modes
Nn = Nm;
Ndam = 12;           % 12 Damage zones

% 12 Mode Shapes Definition
phi = [sin(1*pi*x/a)*sin(1*pi*y/b), sin(1*pi*x/a)*sin(2*pi*y/b), ...
       sin(2*pi*x/a)*sin(1*pi*y/b), sin(2*pi*x/a)*sin(2*pi*y/b), ...
       sin(1*pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(1*pi*y/b)];
       % sin(2*pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(2*pi*y/b), ...
       % sin(3*pi*x/a)*sin(3*pi*y/b), sin(1*pi*x/a)*sin(4*pi*y/b), ...
       % sin(4*pi*x/a)*sin(1*pi*y/b), sin(2*pi*x/a)*sin(4*pi*y/b)];

dxxphi = diff(phi, x, 2);
dxyphi = diff(diff(phi, x), y);
dyyphi = diff(phi, y, 2);

% 12 Damage Zone Coordinates (4x3 Grid over the plate)
damdox = [0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a; ...
          0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a; ...
          0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a];

damdoy = [0 b/3; 0 b/3; 0 b/3; 0 b/3; ...
          b/3 2*b/3; b/3 2*b/3; b/3 2*b/3; b/3 2*b/3; ...
          2*b/3 b; 2*b/3 b; 2*b/3 b; 2*b/3 b];

% 4 Sensor Coordinates
% %             z1      z2          z3      z4        z5    z6       z7       z8       z9      z10     z11      z12     BM      MM        TM
% sensox = [  a/8      3*a/8     5*a/8    7*a/8     a/8    3*a/8    5*a/8    7*a/8    a/8             5*a/8                     a/2       a/2        ];%    a/2                 ];
% sensoy = [  b/6       b/6        b/6      b/6     b/2      b/2     b/2       b/2    5*b/6            5*b/6                     b/2       5*b/6      ];%   b/6            ];
% 

%% sensor from the highest observibility 
sensox = [ 0.1220   0.4780   0.1220  0.4780  0.4881  0.1119  0.2949  0.3051  0.3966  0.2034   0.3254  0.2136 ];
sensoy = [ 0.0718  0.0718  0.3282   0.3282  0.2359  0.2359  0.3282   0.0718  0.1538   0.1538   0.2359  0.2564 ];

%% 3. Matrix Pre-allocation & Symbolic Integration
Mgm = zeros(Nm, Nm);
Kgk = zeros(Nm, Nm);
Fgf = zeros(Nm, 1);
Kd = zeros(Ndam, Nm, Nm);

for i = 1:Nm
    for j = 1:Nm
        Mgm(i,j) = double(ro * h * int(int(phi(i)*phi(j), x, 0, a), y, 0, b));
        Kgk(i,j) = double(D * int(int(dxxphi(i)*dxxphi(j) + 2*dxyphi(i)*dxyphi(j) + dyyphi(i)*dyyphi(j), x, 0, a), y, 0, b));
    end
    Fgf(i,1) = double(int(int(p*phi(i), x, 0, a), y, 0, b));
end

Cgc = 0.0003 * Mgm + 0.0003 * Kgk;     

for i = 1:Ndam
    for j = 1:Nm
        for k = 1:Nm
            Kd(i,j,k) = double(D * int(int(dxxphi(j)*dxxphi(k) + 2*dxyphi(j)*dxyphi(k) + dyyphi(j)*dyyphi(k), x, damdox(i,1), damdox(i,2)), y, damdoy(i,1), damdoy(i,2)));
        end
    end
end

%% 4. Initial Conditions & Damage Scenarios
gamma = 1/2;
beta = 1/2;
x_st = zeros(Nm,1);
v_st = zeros(Nm,1);
acc = zeros(Nm,1);
ac = zeros(Nm,1);
Fg = zeros(Nm,1);
Kgun = Kgk;  

% Damage Injection (Spanning up to 12 zones)
damf = zeros(Ndam, N);
for i = 200:N,  damf(1,i) = (i-200)^.7 * 0.3 / (N-200)^.7; end
for i = 1200:N, damf(2,i) = (i-1200) * 0.6 / (N-1200); end
for i = 600:N,  damf(3,i) = (i-600) * 0.3 / (N-600); end
% for i = 800:N,  damf(4,i) = (i-800) * 0.7 / (N-800); end
% for i = 900:N,  damf(7,i) = (i-900) * 0.4 / (N-900); end
for i = 30:N,   damf(5,i) = (i-30) * 0.2 / (N-30); end
for i = 800:N,  damf(8,i) = (i-800) * 0.4 / (N-800); end
for i = 1500:N, damf(11,i) = (i-1500) * 0.25 / (N-1500); end

%% 5. Sensor Matrix Configuration (12 Modes mapped)
dim = 2 * Nm;
Phi = zeros(nObserv, Nm);
for i = 1:nObserv
    x1 = sensox(i);
    y1 = sensoy(i);
    Phi(i,:) = [sin(1*pi*x1/a)*sin(1*pi*y1/b), sin(1*pi*x1/a)*sin(2*pi*y1/b), ...
                sin(2*pi*x1/a)*sin(1*pi*y1/b), sin(2*pi*x1/a)*sin(2*pi*y1/b), ...
                sin(1*pi*x1/a)*sin(3*pi*y1/b), sin(3*pi*x1/a)*sin(1*pi*y1/b), ...
               ];
end

X = zeros(nObserv, Nm);
%Phi = eye(nObserv);
H = [Phi X ; X Phi];
Haf = Phi;

%% 6. Filter Tuning Matrices (Generalized for 12 Modes/Zones)
% Note: Multipliers retained from your code, structured for 12 dimensions
P = 1e-10 * eye(2*Nm);   
Pdkf2 = 1e-10 * eye(2*Nm); 

Pf = 1e-8 * blkdiag(0.1246, 0.0655,0.0535, 0.1486,0.1804,0.0217, 0.0227,0.1692,0.1699,0.081, 0.0808, 0.1641); 
Pf2 = 1e-8 * blkdiag(0.1246, 0.0655,0.0535, 0.1486,0.1804,0.0217, 0.0227,0.1692,0.1699,0.081, 0.0808, 0.1641);

Q = eye(dim) / 1e-10;  
Qdkf2 = eye(dim) / 1e-10; 

%Pekf = 1e19 * eye(2*Nm + Ndam); 
Pekf = 1e5*blkdiag(1.3421e-14,1.2625e-16, 7.7437e-16, 4.8198e-17,1.0609e-17,7.7538e-17,5.5859e-09, 4.6491e-11,2.8968e-10 ,1.7384e-11, 5.2969e-12, 2.8044e-11, 1.9752e-08,6.7107e-09,6.8089e-09,2.0973e-08,2.2424e-08,3.8958e-09,3.9001e-09, 2.3562e-08,1.4734e-08,9.9054e-09 ,9.9128e-09,1.4875e-08 );

Qekf = 1e-8 * eye(2*Nm + Ndam);
Qf = Pf; 
Qf2 = Pf2; 

R = 10^-1 * eye(2*nObserv);  
Rf = 10^-0.2 * eye(nObserv);
Rdkf2 = 10^-1 * eye(2*nObserv); 
Rf2 = 10^-0.2 * eye(nObserv);
Rekf = 10^-2 * eye(nObserv);

%% 7. Variable Initialization for Time Loop
xx = zeros(dim, 1);
Fo(1) = 0;
FF = zeros(Nm, 1);
aa = zeros(Nm, 1);
vv = zeros(Nm, 1);
fac = zeros(Ndam, 1);
facc = fac;
fd_dkf2 = zeros(Ndam, 1);
fd_est_ex = zeros(Ndam, N);
fd_est_dkf2 = zeros(Ndam, N);
xxp = zeros(2*Nm, 1);
xxp2 = zeros(2*Nm, 1);
xxp_exfd = zeros(2*Nm+Ndam, 1);  
vvp = zeros(Nm, 1);
aap = zeros(Nm, 1);
accp = zeros(Nm, 1);
acp = zeros(Nm, 1);
acp_2 = zeros(Nm, 1);
xp = zeros(Nm, 1);
vp = zeros(Nm, 1);
xrip = zeros(Nm, 1);
vrip = zeros(Nm, 1);
accpri = zeros(Nm, 1);
xx_store_true = zeros(Nm,N);
xx_store_dkf1 = zeros(Nm,N);
xx_store_dkf2 = zeros(Nm, N);
xx_store_ekf = zeros(Nm, N);

Obsv_state_s1 = zeros(2*Nm, 1);
Obsv_fd_s2 = zeros(Ndam, 1);
Obsv_fd_s1 = zeros(Ndam, 1);
Obsv_state_s2 = zeros(2*Nm, 1);
Obsv_ekf = zeros(2*Nm + Ndam, 1);

tem = zeros(Nm, Nm);
xi = zeros(Nm, Ndam);
damfpredict = zeros(N, Ndam);

%% 8. Main Time Integration Loop
for i = 2:N
    % True System Update
    Kgk = Kgun;
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        Kgk = Kgk - damf(j,i) * tem;
    end

    AA = Mgm + Cgc*deltat*gamma + Kgk*beta*deltat^2;
    BB = Cgc*deltat*(1-gamma) + Kgk*(1-2*beta)*deltat^2/2;
    CC = Cgc + deltat*Kgk;

    t = (i-1)*deltat;
    Fg(:) = Fgf * sin(w_1*t);
    % Using \ instead of inv() for numerical stability
    acc(:) = AA \ (Fg(:) - BB*accp(:) - CC*v_st(:) - Kgk*x_st(:));
    accp = acc;

    for kk = 1:Nm
        ac(kk) = acc(kk) + 1e-4*(randn-.5)/50;
    end
    accri = ac + 1e-4*(randn-.5)/70;
    accri(1) = ac(1) + 1e-4*(rand-0.5)/5;

    v_st(:) = vp + deltat*(1-gamma)*accpri + deltat*gamma*accri;
    x_st(:) = xp + deltat*vp + (1-2*beta)/2*deltat^2*accpri + beta*deltat^2*accri;
    vp = v_st;
    xp = x_st;

    for j = 1:Nm
        vri(j) = vrip(j) + (acp(j)+ac(j))/2 * deltat;
        xri(j) = xrip(j) + (vrip(j)+vri(j))/2 * deltat;
    end

    xx_store_true(:,i)=xri;

    z = [(Phi*(xri+xri*0.02*rand)'); (Phi*(vri+vri*0.02*rand)')]';
    zf =( Phi * (ac + ac*0.02*rand))';

    %% Filter Prediction & Update
    disp = xxp(1:Nm) + xxp(Nm+1:2*Nm)*deltat;
    vel = xxp(Nm+1:2*Nm) + acp(1:Nm)*deltat;
    fac = facc;
    Kdx = zeros(Nm, Ndam);

    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        xi(:,j) = tem * disp;
        Kdx(:,j) = xi(:,j);
    end

    F = Fgf * sin(w_1*t);    
    F1 = F;

    Af = eye(Ndam);
    Hf = Haf * (Mgm \ Kdx);
    Cgcx = Mgm \ (Kgun*disp + Cgc*vel);
    HBf = Haf * ((Mgm \ F) - Cgcx);
    O_s1_fd = obsv(Af, Hf);

    % Scheme 2 Stiffness
    Kgks2 = Kgun;
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        Kgks2 = Kgks2 - fd_dkf2(j) * tem; 
    end

    Kdx2 = zeros(Nm, Ndam);
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        xi(:,j) = tem * (xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2);
        Kdx2(:,j) = xi(:,j);
    end

    InvMatS2 = Mgm + deltat*0.5*Cgc + deltat^2*0.25*Kgks2;
    Hdkf2 = Phi * (InvMatS2 \ Kdx2);
    a_dkf2 = Hdkf2 * fd_dkf2 + Phi * (InvMatS2 \ (F1 - Kgun*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2) - Cgc*(xxp2(Nm+1:2*Nm) + deltat^2*0.5*acp_2)));
    O_s2_fd = obsv(Af, Hdkf2);

    % Update Damage Factors (Scheme 1)
    Pnf = Af * Pf * Af' + Qf;
    Kf = Pnf * Hf' / ((Hf * Pnf * Hf') + Rf);
    fac = facc + Kf * (zf' - Hf*fac(:) - HBf);

    % Update Damage Factors (Scheme 2)
    Pnf2 = Af * Pf2 * Af' + Qf2;
    Kdkf2 = Pnf2 * Hdkf2' / ((Hdkf2 * Pnf2 * Hdkf2') + Rf2);
    fd_dkf2 = fd_dkf2 + Kdkf2 * (zf' - a_dkf2);

    %% 2nd Loop for DEKF
    Kgk1 = Kgun;
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        Kgk1 = Kgk1 - fac(j)*tem;  
    end

    Ac = [zeros(Nm,Nm) eye(Nm,Nm); -(Mgm \ Kgk1) -(Mgm \ Cgc)];
    Bc = [zeros(Nm,Nm); Mgm \ eye(Nm)];

    % Scheme 2 Matrix
    Kgdkf2 = Kgun;
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        Kgdkf2 = Kgdkf2 - fd_dkf2(j)*tem;  
    end

    Adkf2 = [zeros(Nm,Nm) eye(Nm,Nm); -(Mgm \ Kgdkf2) -(Mgm \ Cgc)];
    Bdkf2 = [zeros(Nm,Nm); Mgm \ eye(Nm)];

    % Extended Kalman Filter Matrix
    S_ekf = zeros(Nm, Ndam);
    for j = 1:Ndam
       for k = 1:Nm
           for m = 1:Nm
               tem(k,m) = Kd(j,k,m);
           end
       end
       S_ekf(:,j) = tem * xxp_exfd(1:Nm);
    end

    Aex = [zeros(Nm,Nm) , eye(Nm,Nm) , zeros(Nm,Ndam);
           -(Mgm \ Kgk1) , -(Mgm \ Cgc) , (Mgm \ S_ekf);
           zeros(Ndam,Nm), zeros(Ndam,Nm), eye(Ndam,Ndam)];
    Bex = [zeros(Nm,Nm); Mgm \ eye(Nm); zeros(Ndam,Nm)];

    % Discretization Scheme 1
    A = zeros(2*Nm, 2*Nm);
    BB = zeros(2*Nm, 2*Nm);
    for jj = 0:10
        A = A + (Ac*deltat)^jj / factorial(jj);
        BB = BB + (Ac*deltat)^jj / factorial(jj+1);
    end
    B = BB * Bc * deltat;

    % Discretization Scheme 2
    Ax2 = zeros(2*Nm, 2*Nm);
    Bx2_ = zeros(2*Nm, 2*Nm);
    for jj = 0:10
        Ax2 = Ax2 + (Adkf2*deltat)^jj / factorial(jj);
        Bx2_ = Bx2_ + (Adkf2*deltat)^jj / factorial(jj+1);
    end
    Bx2 = Bx2_ * Bdkf2 * deltat;

    % Discretization EKF
    Ax = zeros(2*Nm+Ndam, 2*Nm+Ndam);
    BBx = zeros(2*Nm+Ndam, 2*Nm+Ndam);
    for jj = 0:10
        Ax = Ax + (Aex*deltat)^jj / factorial(jj);
        BBx = BBx + (Aex*deltat)^jj / factorial(jj+1);
    end
    Bx = BBx * Bex * deltat;

    % EKF Update
    xx_exfd = Ax * xxp_exfd + Bx * F;
    Pekfn = Ax * Pekf * Ax' + Qekf; 

    Kgk2 = Kgun;
    for j = 1:Ndam
        for k = 1:Nm
            for m = 1:Nm
                tem(k,m) = Kd(j,k,m);
            end
        end
        Kgk2 = Kgk2 - xx_exfd(2*Nm+j) * tem;  
    end

    for j = 1:Ndam
       for k = 1:Nm
           for m = 1:Nm
               tem(k,m) = Kd(j,k,m);
           end
       end
       S_ekf(:,j) = tem * xx_exfd(1:Nm);
    end

    Hekf = [-Phi*(Mgm \ Kgk2) , -Phi*(Mgm \ Cgc) , Phi*(Mgm \ S_ekf) ];
    z_diff_ekf = zf' - Phi*(Mgm \ (F1 - Kgk2*xx_exfd(1:Nm) - Cgc*xx_exfd(Nm+1 :2*Nm)));
    Kekf = Pekfn * Hekf' / ((Hekf * Pekfn * Hekf') + Rekf);
    xxp_exfd = xx_exfd + Kekf * z_diff_ekf;

    xx_store_ekf(:,i) = xxp_exfd(1:Nm); 
    fd_est_ex(:,i) = xxp_exfd(2*Nm +1: 2*Nm + Ndam);
    O_ex = obsv(Ax, Hekf);

    % Scheme 1 State Update
    xx = A*xxp + B*F;
    Pn = A*P*A' + Q;    
    K = Pn * H' / ((H*Pn*H') + R);
    xx = xx + K*(z' - H*xx);    
    xx_store_dkf1(:,i) = xx(1:Nm);
    vv = (xx(1:Nn) - xxp(1:Nn)) / deltat;
    aa = (vv(1:Nn) - vvp(1:Nn)) / deltat;
    O_s1_st = obsv(Ac, H);

    % Scheme 2 State Update
    xx2 = Ax2*xxp2 + Bx2*F;
    Pndkf2 = Ax2*Pdkf2*Ax2' + Qdkf2;    
    Kx2 = Pndkf2 * H' / ((H*Pndkf2*H') + Rdkf2);
    xx2 = xx2 + Kx2*(z' - H*xx2);   
    xx_store_dkf2(:,i) = xx2(1:Nm);
    fd_est_dkf2(:,i) = fd_dkf2;
    O_s2_st = obsv(Adkf2, H);

    %% Accumulating Observability Metrics
    for z_idx = 1:Ndam
        Obsv_fd_s1(z_idx) = Obsv_fd_s1(z_idx) + norm(O_s1_fd(:, z_idx))^2;
        Obsv_fd_s2(z_idx) = Obsv_fd_s2(z_idx) + norm(O_s2_fd(:, z_idx))^2;
    end
    for k = 1:2*Nm
        Obsv_state_s1(k) = Obsv_state_s1(k) + norm(O_s1_st(:, k))^2;
        Obsv_state_s2(k) = Obsv_state_s2(k) + norm(O_s2_st(:, k))^2;
    end
    for k = 1:(2*Nm+Ndam)
        Obsv_ekf(k) = Obsv_ekf(k) + norm(O_ex(:,k))^2;
    end

    damfpredict(i,:) = fac;
    facc = fac;

    % Covariance Updates
    P = (eye(dim) - K*H) * Pn;
    Pf = (eye(Ndam) - Kf*Hf) * Pnf;
    Pdkf2 = (eye(dim) - Kx2*H) * Pndkf2;
    Pf2 = (eye(Ndam) - Kdkf2*Hdkf2) * Pnf2;

    xxp = xx;
    vvp = vv;
    aap = aa;
    acp_2 = (xx2(Nm+1:2*Nm) - xxp2(Nm+1:2*Nm)) / deltat;
    xxp2 = xx2;
    accp = acc;
    acp = ac;
    xrip = xri;
    vrip = vri;
    accpri = accri;
end

%% 9. Visualization & Plotting (Updated for 12 Damages)
figure(1);
sgtitle('Damage Factors: True vs Estimated (Zones 1-12)');
for i = 1:Ndam
    subplot(4,3,i)

    plot(damf(i,:), 'k-'); hold on;
    plot(damfpredict(:,i), 'b--');
    plot(fd_est_dkf2(i,:), 'm--');
    plot(fd_est_ex(i,:), 'g--');
    title(sprintf('Zone %d', i));
    xlabel('time');
    ylabel('Damage factor');
    hold off
end

%% 10. Secondary Observability Analysis
dt = 0.00001;       
T_obs = 0.1;    
T = round(T_obs/dt);
q = zeros(Nm,N);
qd = zeros(Nm,N);

for j = 1:Ndam
    for i = 1:Nm
        for k = 1:Nm
            KKd(i,k) = Kd(j,i,k);
        end
    end
    K_zone{j} = KKd;
end

F = Fgf * sin(w_1*(0:N-1)*deltat);

for n = 2:N
    qdd = Mgm \ (F(:,n-1) - Kgk*q(:,n-1));
    qd(:,n) = qd(:,n-1) + deltat*qdd;
    q(:,n) = q(:,n-1) + deltat*qd(:,n-1);
end

for i = 1:Ndam
    G = eye(Nm)' * K_zone{i} * eye(Nm);      
    H_at = Phi * (Mgm \ G);         
    for n = 1:N
        D_i(:,n) = H_at * q(:,n);       
    end
    O(i) = mean(vecnorm(D_i,2,1));      
end

for i = 1:2*Nm
    Hs = [Phi , X ; X , Phi];         
    for n = 1:N
        DS_i(:,n) = Hs * [q(:,n); qd(:,n)];     
    end
    OS(i) = mean(vecnorm(DS_i,2,1));      
end

OS = OS / max(OS);          
O = O / max(O);          

%% 11. Final Plotting (Expanded for 12/24 ranges)
figure(2);
bar(OS);
xlabel('No. of modes (disp & Vel: 1-24)');
ylabel('Normalized observability');
title('State Observability');
grid on;

figure(3);
bar(O);
xlabel('Damage zone (1-12)');
ylabel('Normalized observability');
title('Zone-wise damage observability');
grid on;

% Label Generics for 12/24 sizes
zone_labels = arrayfun(@(x) sprintf('Z%d', x), 1:12, 'UniformOutput', false);
mode_labels = arrayfun(@(x) sprintf('m%d', x), 1:24, 'UniformOutput', false);

figure(4)
bar(Obsv_fd_s1);
set(gca, 'XTickLabel', zone_labels);
title('Observability (Fd) over Time for scheme 1');

figure(5)
bar(Obsv_fd_s2);
set(gca, 'XTickLabel', zone_labels);
title('Observability (Fd) over Time for scheme 2');

figure(6)
bar(Obsv_state_s1);
set(gca, 'XTickLabel', mode_labels);
title('Observability Strength (modes) over Time for scheme 1');

figure(7)
bar(Obsv_state_s2);
set(gca, 'XTickLabel', mode_labels);
title('Observability Strength (modes) over Time for scheme 2');

figure(8)
bar(Obsv_ekf(2*Nm+1:2*Nm+Ndam));
set(gca, 'XTickLabel', zone_labels);
title('Observability Strength (dams) over Time for EKF');

P1_s1 = 1./Obsv_state_s1;
P1_s2 = 1./Obsv_state_s2;

P2_s1 = 1./Obsv_fd_s1;
P2_s2 = 1./Obsv_fd_s2;
%%
figure(10);
  for i =1:Nm
        subplot(2,3,i);
        plot(linspace(0,deltat*N,N),xx_store_true(i,:),'k-'); hold on;
        plot(linspace(0,deltat*N,N),xx_store_dkf2(i,:),'m--'); hold on; 
        plot(linspace(0,deltat*N,N),xx_store_dkf1(i,:),'b--'); hold on;
        plot(linspace(0,deltat*N,N),xx_store_ekf(i,:),'g--');
        xlabel('time');
        ylabel('displacement');
        title('displacement vs time');

  end


%% 
% P_ekf = 1./Obsv_ekf;
% for i=1:24
%     P_ekf(i)
% end

%% cHECKING THE REDUNDANCY 
% [R, pivots] = rref(Phi);
% 
% % 'pivots' will list the columns that are independent.
% % Any number from 1 to 12 NOT in 'pivots' is a redundant column.
% redundant_cols = setdiff(1:12, R);
% disp('Redundant Columns are:');
% disp(redundant_cols);


