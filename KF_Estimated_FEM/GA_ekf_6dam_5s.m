% % %% Initial 
%
%     first run normal for the systen save imp data in physics plate etc then fetch them for the rest of 
%      GA algorithm 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % clear all; close all;
% % syms tau
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
% %     Bx=zeros(2*Nm+Ndam,Nm);
% %     for jj=0:10
% %         Ax=Ax+(Aex*deltat)^jj/factorial(jj);
% %         BBx=BBx+(Aex*deltat)^jj/factorial(jj+1);
% %     end
% % 
% %     Bx=BBx*Bex*deltat;
% %      xx_exfd = Ax*xxp_exfd + Bx*F;
% %     Pekfn = Ax*Pekf*transpose(Ax)+ Qekf; 
% %      Kgk2 = Kgun;
% % 
% %      for j=1:Ndam
% %         for k=1:Nm
% %             for m=1:Nm
% %                 tem(k,m)=Kd(j,k,m);
% %             end
% %         end
% %         Kgk2=Kgk2-xx_exfd(2*Nm+j)*tem;  % here must be fac or facc damf will not be there
% %     end
% %       S_ekf=zeros(Nm,Ndam);
% %     for j=1:Ndam
% %        for k=1:Nm
% %            for m=1:Nm
% %                tem(k,m)=Kd(j,k,m);
% %            end
% %        end
% %        S_ekf(:,j) = tem*xx_exfd(1:Nm);
% %     end
% % 
% % 
% % 
% %     Hekf = [-Phi*inv(Mgm)*Kgk2 , -Phi*inv(Mgm)*Cgc , Phi*inv(Mgm)*S_ekf ];
% %     z_diff_ekf = zf' - Phi*( Mgm\(F1 - Kgk2*xx_exfd(1:Nm) - Cgc*xx_exfd(Nm+1 :2*Nm) )) ;
% %     Kekf = Pekfn*transpose(Hekf)*inv((Hekf*Pekfn*transpose(Hekf)) + Rekf);
% %     xxp_exfd = xx_exfd + Kekf*z_diff_ekf;
% % 
% %     % Storing the Values of the updated ekf 
% %     xx_store_ekf = xxp_exfd(1:2*Nm); %(  first Nm 'u' then Nm 'vel')
% %     fd_est_ex(:,i) = xxp_exfd(2*Nm +1: 2*Nm + Ndam);
% % 
% % 
% %    % 0*B*Fg(:,i-1)
% %    % For Scheme 1 
% %     xx=A*xxp+B*F;
% %     Pn=A*P*transpose(A)+Q;    
% %     K=Pn*transpose(H)*inv((H*Pn*transpose(H))+R);
% %     xx=xx+K*(z'-H*xx);    
% %     vv=(xx(1:Nn)-xxp(1:Nn))/deltat;
% %     aa=(vv(1:Nn)-vvp(1:Nn))/deltat;
% %     FF=Mgm*aa(:)+Cgc*vv(:)+Kgk*xx(1:Nn);
% %     %KK(:,i)=K;
% % 
% % 
% % 
% %     %Fo(i)=Fo(i-1)+deltat/2*(xx(5,i-1)+xx(5,i));
% % 
% %     % For Scheme 2 
% %     xx2=Ax2*xxp2+Bx2*F;
% %     Pndkf2=Ax2*Pdkf2*transpose(Ax2)+Qdkf2;    
% %     Kx2=Pndkf2*transpose(H)*inv((H*Pndkf2*transpose(H))+Rdkf2);
% %     xx2=xx2+Kx2*(z'-H*xx2);   
% %     xx_store_dkf2(:,i) = xx2;
% %     fd_est_dkf2(:,i) = fd_dkf2;
% %     % vv2=(xx2(1:Nn)-xxp2(1:Nn))/deltat;
% % 
% % 
% %  reba=1000;
% % 
% %    % fprintf(Fil1,'%f\n',fac');
% % 
% %       damfpredict(i,:)=fac;
% %    facc=fac;
% % %     
% % %     end
% %     P=(eye(dim)-K*H)*Pn;
% %     Pf=(eye(Ndam)-Kf*Hf)*Pnf;
% % 
% %      Pdkf2=(eye(dim)-Kx2*H)*Pndkf2;
% %     Pf2=(eye(Ndam)-Kdkf2*Hdkf2)*Pnf2;
% % 
% % 
% %     xxp=xx;
% %     vvp=vv;
% %     aap=aa;
% % 
% %     acp_2 = (xx2(Nm+1:2*Nm)-xxp2(Nm+1:2*Nm))/deltat;
% %     xxp2 = xx2;
% % 
% % 
% %     accp=acc;
% %     acp=ac;
% %     xp=x;
% %     vp=v;
% %     xrip=xri;
% %     vrip=vri;
% %     accpri=accri;
% % 
% % end
% % 
% % %% 
% % for i=1:6
% %     subplot(2,3,i)
% %     plot(damf(i,:),'k-');
% %     hold on;
% %     plot(damfpredict(:,i),'b--'); hold on;
% %     plot(fd_est_dkf2(i,:),'m--');
% %     plot(fd_est_ex(i,:),'g--');
% %     hold off
% % end
% % 
% 
% 
% 
% 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %
% % %    %6 Modes, 12 Damage Factors 8SENSORS - Dual KF & EKF with Observability
% % %
% % %
% % %%%%% _____________________________________________________________________________________________________________________________________________________________________________________________
% % 
% clear all;
% close all; clc;
% p=1000;
% ro=2700; % ( in kg/m3)
% yo =70*10^9;%youngs modulus
% neu=0.3;
% h=0.001;%thickness
% a=0.6; %( m)
% b=0.4; % (m)
% divx = 3;  %input(' enter the no of divisions on the breadth ');
% % ex=length of element in x direction
% ex =a/(divx);
% divy=3;  %input(' enter the no of divisions on the length ');
% % ey=length of element in y direction
% ey =b/(divy);
% %% Natural frequencies for the above 
% %% NB: when frequencies are around 8,9 w1 we find scheme 2 is better than the other model
% w_1 = 100*pi;   %1*137.2892;
% w_6 =  897.6603; 
% 
% N = 50000;
% deltat = 0.00005;
% nObserv= 5 ;       % No of observable zones ( sensor is acting) 
% 
% syms x y
% 
% 
% 
% Nm=6;
% Nn=Nm;
% %phi=[sin(pi*x/a)*sin(pi*y/b) sin(3*pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(3*pi*y/b) sin(3*pi*x/a)*sin(3*pi*y/b) sin(4*pi*x/a)*sin(4*pi*y/b)];
% phi=[sin(pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(2*pi*y/b) sin(2*pi*x/a)*sin(pi*y/b) sin(2*pi*x/a)*sin(2*pi*y/b) sin(1*pi*x/a)*sin(3*pi*y/b) sin(3*pi*x/a)*sin(1*pi*y/b)];
% %phi=[sin(pi*x/a)*sin(pi*y/b)  sin(3*pi*x/a)*sin(pi*y/b) sin(pi*x/a)*sin(3*pi*y/b)  sin(3*pi*x/a)*sin(3*pi*y/b)  sin(pi*x/a)*sin(5*pi*y/b)  sin(5*pi*x/a)*sin(1*pi*y/b)];
% 
% dxxphi=diff(phi,x,2);
% dxyphi=diff(diff(phi,x),y);
% dyyphi=diff(phi,y,2);
% D=(yo*h^3)/(12*(1-neu^2));
% 
% Ndam=6;  % Total no. of damage zones ( each having own fd) not the sensors 
% % damdox=[0 a/4;a/4 a/2; a/2 3*a/4;3*a/4 a;0 a/4;a/4 a/2];
% % damdoy=[0 b/4;0 b/4;0 b/4;0 b/4;b/4 b/2;b/4 b/2;b/4 b/2];
% 
% damdox=[0 a/3; a/3 2*a/3;2*a/3 a;0 a/3; a/3 2*a/3;2*a/3 a];
% damdoy=[0 b/2;0 b/2;0 b/2;b/2 b;b/2 b;b/2 b];
% 
% % Most obsv Sensor location for 6 sensor placement 
% % sensox = [0.2949  0.1627  0.4271  0.1525  0.4475  0.3051];
% % sensoy = [ 0.1949 0.1026  0.1026  0.2769  0.2769  0.3179];
% 
% 
% 
% 
% % Middle term is added as the 
%  % sensox=[a/6 a/2 5*a/6  a/6  5*a/6   a/2 ];
%  % sensoy=[b/4 b/4 b/4  3*b/4  3*b/4  b/2 ];
% 
% 
% % %  % Effective for 5 sensor location 
% %     sensox = [0.2949  0.1627  0.4271  0.1525   ];
% %     sensoy = [ 0.1949 0.1026  0.1026  0.2769  ];
%  sensox=[  a/6   5*a/6   1*a/6   5*a/6   a/2  ];
%  sensoy=[   b/4   b/4    3*b/4    3*b/4   b/2 ];
% 
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
% Kgun=Kgk;  % kGK Is the undamaged one 
% Haf=[1 0 0 0];
% 
% damf=zeros(Ndam,N);
% 
%   for i=200:N
%       damf(1,i)=(i-200)^.7*0.3/(N-200)^.7;
%   end
% % % 
%   for i=1200:N
%       damf(2,i)=(i-1200)*0.6/(N-1200);
%   end
%  
%   for i=600:N
%      damf(3,i)=(i-600)*0.12/(N-600);
%  end
% 
% % for i=10:N
% %     damf(4,i)=((i-10)^0.99)*0.3/(N-10);
% % end
% % 
% for i=30:N
%     damf(5,i)=(i-30)*0.45/(N-30);
% end
% 
% 
% for i=300:N
%     damf(6,i)=(i-300)*0.2/(N-300);
% end
% 
% % % % NON LINEAR DAMAGES 
% %   for i=2000:N
% %       damf(1,i)=(i-2000)^.7*0.3/(N-2000)^.7;
% %   end
% % % % 
% %   for i=12000:N
% %       damf(2,i)=(i-12000)*0.6/(N-12000);
% %   end
% % %  
% %   for i=600:N
% %      damf(3,i)=((i-600)^0.8)*0.12/(N-600);
% %  end
% % 
% % for i=100:N
% %     damf(4,i)=(i-100)*0.3/(N-100);
% % end
% % 
% % for i=300:N
% %     damf(5,i)=(i-300)*0.2/(N-300);
% % end
% % 
% % 
% % for i=300:N
% %     damf(6,i)=(i-300)*0.2/(N-300);
% % end
% 
% 
% 
% dim=2*Nm;
% 
% 
% 
% for i=1:nObserv
%     x1=sensox(i);
%     y1=sensoy(i);
%    % Phi(i,:)=[sin(pi*x1/a)*sin(pi*y1/b) sin(2*pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(2*pi*y1/b) sin(2*pi*x1/a)*sin(2*pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b) sin(3*pi*x1/a)*sin(3*pi*y1/b)  sin(pi*x1/a)*sin(4*pi*y1/b) sin(2*pi*x1/a)*sin(4*pi*y1/b) sin(3*pi*x1/a)*sin(4*pi*y1/b) sin(4*pi*x1/a)*sin(1*pi*y1/b) sin(4*pi*x1/a)*sin(2*pi*y1/b) sin(4*pi*x1/a)*sin(3*pi*y1/b)];
%      %Phi(i,:) = [sin(pi*x1/a)*sin(pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b)  sin(3*pi*x1/a)*sin(3*pi*y1/b)  sin(pi*x1/a)*sin(5*pi*y1/b)  sin(5*pi*x1/a)*sin(1*pi*y1/b)];
%      Phi(i,:) = [ sin(pi*x1/a)*sin(pi*y1/b) sin(1*pi*x1/a)*sin(2*pi*y1/b) sin(2*pi*x1/a)*sin(1*pi*y1/b) sin(2*pi*x1/a)*sin(2*pi*y1/b) sin(pi*x1/a)*sin(3*pi*y1/b) sin(3*pi*x1/a)*sin(pi*y1/b)];
% end
% maxVal = max(Phi, [], 'all');
% Phi = Phi./maxVal;
% 
% X=zeros(nObserv,Nm);
% 
% H=[Phi  X ; X  Phi];
% Haf=Phi;
% 
% %% 6. Filter Tuning Matrices (Generalized for 12 Modes/Zones)
% % Note: Multipliers retained from your code, structured for 12 dimensions
% P=[ones(Nm,1)*10^(-4);ones(Nm,1)*10^(-4)];   % 10^-10  for in 4 sensors
%  P=diag(P);
%  Pdkf2=P;
% %P = 1e-5*((3.8302e-61)^-1)*blkdiag(1.3013e-65,2.6286e-69, 3.4733e-66, 2.8234e-70, 1.6619e-73,6.9362e-69 , 3.2922e-61 , 6.7425e-66, 3.8302e-61, 3.1583e-65, 5.8089e-70, 1.1842e-63);   % for 6 sensor and usual mode shapes
% %Pdkf2= 1e-5*((3.8302e-61)^-1)*blkdiag( 1.3047e-65,2.6293e-69, 3.4802e-66, 2.8268e-70, 1.6648e-73,6.9442e-69 , 3.2971e-61 , 6.7527e-66, 3.8364e-61, 3.1618e-65, 5.8148e-70, 1.1854e-63); % for 6 sensor and usual mode shapes
% 
% %Pf=eye(Ndam)*10000;  % BEST AT 10000 
% Pf = 1e-12*blkdiag(0.3728,0.3122,0.4581,0.3710,0.4272,0.3886);% for 6 sensor scheme 1;
% Pf2 = 1e-12*blkdiag(0.3739 ,0.3122,0.7081, 0.6620,  0.6101 , 0.7058); % for 6 sensor scheme 2
%  Q=eye(dim)/10^-14;  %  % scheme 1 -> eye(dim)/10^-10 for 6 sensor (unchange modes)
%  Qdkf2=eye(dim)/10^-14; %
% dimf=1;
% 
% %Pekf = 1e19 * eye(2*Nm + Ndam); 
% Q_ekf_fd = 1e-1*blkdiag(5.2941e-11, 8.4701e-11, 1.9810e-11,2.6663e-11 ,2.6273e-11,8.6801e-10 );      %    0.2663 , 0.2135 , 0.3137, 0.1538 ,0.4348, 0.8606);
% P_state_ekf = 1e-2*blkdiag( 1.3499e-15 , 1.6672e-15 , 9.6627e-17 , 4.6831e-18 ,1.2994e-18 , 1.0910e-18,1.1649e-17 ,1.4995e-08, 1.8522e-10, 1.0734e-09,3.6870e-11 , 1.2121e-11 );
% Pekf =  blkdiag(1e0*P_state_ekf,Q_ekf_fd);
% Qekf = 1e-8 * eye(2*Nm + Ndam);
% Qf = Pf; %1e-12*blkdiag( 0.4793  , 0.3464 ,0.5248, 0.4994 ,0.3564, 0.5304);
% Qf2 =Pf2; %    % 1e1*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);
% 
% R=1*eye(2*nObserv)*10^(4);  % for state ( disp & vel) of Scheme - 1 for 6 sen ->1*eye(2*nObserv)*10^(-4) works better ( no noise)  
% Rf=1*eye(nObserv)*10^(2);
% Rf2=1*eye(nObserv)*10^(2);
% %Rf = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);
% 
% Rdkf2=1*eye(2*nObserv)*10^(4); % for state ( disp & vel) of Scheme - 2
% %Rf2 = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);
% 
% Rekf = 10^0 * eye(nObserv);
% 
% %% 7. Variable Initialization for Time Loop
% xx = zeros(dim, 1);
% Fo(1) = 0;
% FF = zeros(Nm, 1);
% aa = zeros(Nm, 1);
% vv = zeros(Nm, 1);
% fac = zeros(Ndam, 1);
% facc = fac;
% fd_dkf2 = zeros(Ndam, 1);
% fd_est_ex = zeros(Ndam, N);
% fd_est_dkf2 = zeros(Ndam, N);
% xxp = zeros(2*Nm, 1);
% xxp2 = zeros(2*Nm, 1);
% xxp_exfd = zeros(2*Nm+Ndam, 1);  
% vvp = zeros(Nm, 1);
% aap = zeros(Nm, 1);
% accp = zeros(Nm, 1);
% acp = zeros(Nm, 1);
% acp_2 = zeros(Nm, 1);
% xp = zeros(Nm, 1);
% vp = zeros(Nm, 1);
% xrip = zeros(Nm, 1);
% vrip = zeros(Nm, 1);
% accpri = zeros(Nm, 1);
% xx_store_true = zeros(Nm,N);
% xx_store_dkf1 = zeros(Nm,N);
% xx_store_dkf2 = zeros(Nm, N);
% xx_store_ekf = zeros(Nm, N);
% x_st = zeros(Nm,1);
% v_st = zeros(Nm,1);
% acc = zeros(Nm,1);
% ac = zeros(Nm,1);
% Fg = zeros(Nm,1);
% Obsv_state_s1 = zeros(2*Nm, 1);
% Obsv_fd_s2 = zeros(Ndam, 1);
% Obsv_fd_s1 = zeros(Ndam, 1);
% Obsv_state_s2 = zeros(2*Nm, 1);
% Obsv_ekf = zeros(2*Nm + Ndam, 1);
% 
% tem = zeros(Nm, Nm);
% xi = zeros(Nm, Ndam);
% damfpredict = zeros(N, Ndam);
% 
% %% sAVING FOR THE GA SOLUTION 
% % Collect everything into a clean structure
% plate_physics.Nm = Nm;
% plate_physics.Ndam = Ndam;
% plate_physics.nObserv = nObserv;
% plate_physics.deltat = deltat;
% plate_physics.w_1 = w_1;
% 
% % Numerical Matrices
% plate_physics.Mgm = double(Mgm);
% plate_physics.Kgun = double(Kgun);
% plate_physics.Cgc = double(Cgc);
% plate_physics.Fgf = double(Fgf);
% plate_physics.Kd = double(Kd); % The 3D array
% plate_physics.Phi = double(Phi);
% 
% % Target Data for GA to match
% plate_physics.true_damage = damf; 
% 
% % Important: Pre-calculate the Inverse of Mgm 
% % Since Mgm is constant, (Mgm \ F) is faster if you pre-calculate invMgm
% plate_physics.invMgm = inv(double(Mgm)); 
% 
% save('PlateData.mat', 'plate_physics');
% disp('Data successfully wrapped for GA optimization.');
% 
% %% 8. Main Time Integration Loop
% for i = 2:N
%     % True System Update
%     Kgk = Kgun;
%     for j = 1:Ndam
%         for k = 1:Nm
%             for m = 1:Nm
%                 tem(k,m) = Kd(j,k,m);
%             end
%         end
%         Kgk = Kgk - damf(j,i) * tem;
%     end
% 
%     AA = Mgm + Cgc*deltat*gamma + Kgk*beta*deltat^2;
%     BB = Cgc*deltat*(1-gamma) + Kgk*(1-2*beta)*deltat^2/2;
%     CC = Cgc + deltat*Kgk;
% 
%     t = (i-1)*deltat;
%     Fg(:) = Fgf * sin(w_1*t);
%     % Using \ instead of inv() for numerical stability
%     acc(:) = AA \ (Fg(:) - BB*accp(:) - CC*v_st(:) - Kgk*x_st(:));
%     accp = acc;
% 
%     for kk = 1:Nm
%         ac(kk) = acc(kk) + 1e-4*(randn-.5)/50;
%     end
%     accri = ac + 1e-4*(randn-.5)/70;
%     accri(1) = ac(1) + 1e-4*(rand-0.5)/5;
% 
%     v_st(:) = vp + deltat*(1-gamma)*accpri + deltat*gamma*accri;
%     x_st(:) = xp + deltat*vp + (1-2*beta)/2*deltat^2*accpri + beta*deltat^2*accri;
%     vp = v_st;
%     xp = x_st;
% 
%     for j = 1:Nm
%         vri(j) = vrip(j) + (acp(j)+ac(j))/2 * deltat;
%         xri(j) = xrip(j) + (vrip(j)+vri(j))/2 * deltat;
%     end
% 
%     xx_store_true(:,i)=xri;
% 
%     z = [(Phi*(xri+xri*0.02*rand)'); (Phi*(vri+vri*0.02*rand)')]';
%     zf =( Phi * (ac + ac*0.03*rand))';
% 
%     %% Filter Prediction & Update
%     disp = xxp(1:Nm) + xxp(Nm+1:2*Nm)*deltat;
%     vel = xxp(Nm+1:2*Nm) + acp(1:Nm)*deltat;
%     fac = facc;
%     Kdx = zeros(Nm, Ndam);
% 
%     for j = 1:Ndam
%         for k = 1:Nm
%             for m = 1:Nm
%                 tem(k,m) = Kd(j,k,m);
%             end
%         end
%         xi(:,j) = tem * disp;
%         Kdx(:,j) = xi(:,j);
%     end
% 
%     F = Fgf * sin(w_1*t);    
%     F1 = F;
% 
%     Af = eye(Ndam);
%     Hf = Haf * (Mgm \ Kdx);
%     Cgcx = Mgm \ (Kgun*disp + Cgc*vel);
%     HBf = Haf * ((Mgm \ F) - Cgcx);
%     O_s1_fd = obsv(Af, Hf);
% 
%     % Scheme 2 Stiffness
%     Kgks2 = Kgun;
%     for j = 1:Ndam
%         for k = 1:Nm
%             for m = 1:Nm
%                 tem(k,m) = Kd(j,k,m);
%             end
%         end
%         Kgks2 = Kgks2 - fd_dkf2(j) * tem; 
%     end
% 
%     Kdx2 = zeros(Nm, Ndam);
%     for j = 1:Ndam
%         for k = 1:Nm
%             for m = 1:Nm
%                 tem(k,m) = Kd(j,k,m);
%             end
%         end
%         xi(:,j) = tem * (xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2);
%         Kdx2(:,j) = xi(:,j);
%     end
% 
%     InvMatS2 = Mgm + deltat*0.5*Cgc + deltat^2*0.25*Kgks2;
%     Hdkf2 = Phi * (InvMatS2 \ Kdx2);
%     a_dkf2 = Hdkf2 * fd_dkf2 + Phi * (InvMatS2 \ (F1 - Kgun*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2) - Cgc*(xxp2(Nm+1:2*Nm) + deltat^2*0.5*acp_2)));
%     O_s2_fd = obsv(Af, Hdkf2);
% 
%     % Update Damage Factors (Scheme 1)
%     Pnf = Af * Pf * Af' + Qf;
%     Kf = Pnf * Hf' / ((Hf * Pnf * Hf') + Rf);
%     fac = facc + Kf * (zf' - Hf*fac(:) - HBf);
% 
%     % Update Damage Factors (Scheme 2)
%     Pnf2 = Af * Pf2 * Af' + Qf2;
%     Kdkf2 = Pnf2 * Hdkf2' / ((Hdkf2 * Pnf2 * Hdkf2') + Rf2);
%     fd_dkf2 = fd_dkf2 + Kdkf2 * (zf' - a_dkf2);
% 
%     %% 2nd Loop for DEKF
%     Kgk1 = Kgun;
%     for j = 1:Ndam
%         for k = 1:Nm
%             for m = 1:Nm
%                 tem(k,m) = Kd(j,k,m);
%             end
%         end
%         Kgk1 = Kgk1 - fac(j)*tem;  
%     end
% 
%     Ac = [zeros(Nm,Nm) eye(Nm,Nm); -(Mgm \ Kgk1) -(Mgm \ Cgc)];
%     Bc = [zeros(Nm,Nm); Mgm \ eye(Nm)];
% 
%     % Scheme 2 Matrix
%     Kgdkf2 = Kgun;
%     for j = 1:Ndam
%         for k = 1:Nm
%             for m = 1:Nm
%                 tem(k,m) = Kd(j,k,m);
%             end
%         end
%         Kgdkf2 = Kgdkf2 - fd_dkf2(j)*tem;  
%     end
% 
%     Adkf2 = [zeros(Nm,Nm) eye(Nm,Nm); -(Mgm \ Kgdkf2) -(Mgm \ Cgc)];
%     Bdkf2 = [zeros(Nm,Nm); Mgm \ eye(Nm)];
% 
%     % Extended Kalman Filter Matrix
%     S_ekf = zeros(Nm, Ndam);
%     for j = 1:Ndam
%        for k = 1:Nm
%            for m = 1:Nm
%                tem(k,m) = Kd(j,k,m);
%            end
%        end
%        S_ekf(:,j) = tem * xxp_exfd(1:Nm);
%     end
% 
%     Aex = [zeros(Nm,Nm) , eye(Nm,Nm) , zeros(Nm,Ndam);
%            -(Mgm \ Kgk1) , -(Mgm \ Cgc) , (Mgm \ S_ekf);
%            zeros(Ndam,Nm), zeros(Ndam,Nm), eye(Ndam,Ndam)];
%     Bex = [zeros(Nm,Nm); Mgm \ eye(Nm); zeros(Ndam,Nm)];
% 
%     % Discretization Scheme 1
%     A = zeros(2*Nm, 2*Nm);
%     BB = zeros(2*Nm, 2*Nm);
%     for jj = 0:10
%         A = A + (Ac*deltat)^jj / factorial(jj);
%         BB = BB + (Ac*deltat)^jj / factorial(jj+1);
%     end
%     B = BB * Bc * deltat;
% 
%     % Discretization Scheme 2
%     Ax2 = zeros(2*Nm, 2*Nm);
%     Bx2_ = zeros(2*Nm, 2*Nm);
%     for jj = 0:10
%         Ax2 = Ax2 + (Adkf2*deltat)^jj / factorial(jj);
%         Bx2_ = Bx2_ + (Adkf2*deltat)^jj / factorial(jj+1);
%     end
%     Bx2 = Bx2_ * Bdkf2 * deltat;
% 
%     % Discretization EKF
%     Ax = zeros(2*Nm+Ndam, 2*Nm+Ndam);
%     BBx = zeros(2*Nm+Ndam, 2*Nm+Ndam);
%     for jj = 0:10
%         Ax = Ax + (Aex*deltat)^jj / factorial(jj);
%         BBx = BBx + (Aex*deltat)^jj / factorial(jj+1);
%     end
%     Bx = BBx * Bex * deltat;
% 
%     % EKF Update
%     xx_exfd = Ax * xxp_exfd + Bx * F;
%     Pekfn = Ax * Pekf * Ax' + Qekf; 
% 
%     Kgk2 = Kgun;
%     for j = 1:Ndam
%         for k = 1:Nm
%             for m = 1:Nm
%                 tem(k,m) = Kd(j,k,m);
%             end
%         end
%         Kgk2 = Kgk2 - xx_exfd(2*Nm+j) * tem;  
%     end
% 
%     for j = 1:Ndam
%        for k = 1:Nm
%            for m = 1:Nm
%                tem(k,m) = Kd(j,k,m);
%            end
%        end
%        S_ekf(:,j) = tem * xx_exfd(1:Nm);
%     end
% 
%     Hekf = [-Phi*(Mgm \ Kgk2) , -Phi*(Mgm \ Cgc) , Phi*(Mgm \ S_ekf) ];
%     z_diff_ekf = zf' - Phi*(Mgm \ (F1 - Kgk2*xx_exfd(1:Nm) - Cgc*xx_exfd(Nm+1 :2*Nm)));
%     Kekf = Pekfn * Hekf' / ((Hekf * Pekfn * Hekf') + Rekf);
%     xxp_exfd = xx_exfd + Kekf * z_diff_ekf;
% 
%     xx_store_ekf(:,i) = xxp_exfd(1:Nm); 
%     fd_est_ex(:,i) = xxp_exfd(2*Nm +1: 2*Nm + Ndam);
%     O_ex = obsv(Ax, Hekf);
% 
%     % Scheme 1 State Update
%     xx = A*xxp + B*F;
%     Pn = A*P*A' + Q;    
%     K = Pn * H' / ((H*Pn*H') + R);
%     xx = xx + K*(z' - H*xx);    
%     xx_store_dkf1(:,i) = xx(1:Nm);
%     vv = (xx(1:Nn) - xxp(1:Nn)) / deltat;
%     aa = (vv(1:Nn) - vvp(1:Nn)) / deltat;
%     O_s1_st = obsv(Ac, H);
% 
%     % Scheme 2 State Update
%     xx2 = Ax2*xxp2 + Bx2*F;
%     Pndkf2 = Ax2*Pdkf2*Ax2' + Qdkf2;    
%     Kx2 = Pndkf2 * H' / ((H*Pndkf2*H') + Rdkf2);
%     xx2 = xx2 + Kx2*(z' - H*xx2);   
%     xx_store_dkf2(:,i) = xx2(1:Nm);
%     fd_est_dkf2(:,i) = fd_dkf2;
%     O_s2_st = obsv(Adkf2, H);
% 
%     %% Accumulating Observability Metrics
%     for z_idx = 1:Ndam
%         Obsv_fd_s1(z_idx) = Obsv_fd_s1(z_idx) + norm(O_s1_fd(:, z_idx))^2;
%         Obsv_fd_s2(z_idx) = Obsv_fd_s2(z_idx) + norm(O_s2_fd(:, z_idx))^2;
%     end
%     for k = 1:2*Nm
%         Obsv_state_s1(k) = Obsv_state_s1(k) + norm(O_s1_st(:, k))^2;
%         Obsv_state_s2(k) = Obsv_state_s2(k) + norm(O_s2_st(:, k))^2;
%     end
%     for k = 1:(2*Nm+Ndam)
%         Obsv_ekf(k) = Obsv_ekf(k) + norm(O_ex(:,k))^2;
%     end
% 
%     damfpredict(i,:) = fac;
%     facc = fac;
% 
%     % Covariance Updates
%     P = (eye(dim) - K*H) * Pn;
%     Pf = (eye(Ndam) - Kf*Hf) * Pnf;
%     Pdkf2 = (eye(dim) - Kx2*H) * Pndkf2;
%     Pf2 = (eye(Ndam) - Kdkf2*Hdkf2) * Pnf2;
% 
%     xxp = xx;
%     vvp = vv;
%     aap = aa;
%     acp_2 = (xx2(Nm+1:2*Nm) - xxp2(Nm+1:2*Nm)) / deltat;
%     xxp2 = xx2;
%     accp = acc;
%     acp = ac;
%     xrip = xri;
%     vrip = vri;
%     accpri = accri;
% end
% 
% %% 9. Visualization & Plotting (Updated for 12 Damages)
% figure(1);
% sgtitle('Damage Factors: True vs Estimated (Zones 1-12)');
% for i = 1:Ndam
%     subplot(4,3,i)
% 
%     plot(damf(i,:), 'k-'); hold on;
%     plot(damfpredict(:,i), 'b--');
%     plot(fd_est_dkf2(i,:), 'm--');
%     plot(fd_est_ex(i,:), 'g--');
%     title(sprintf('Zone %d', i));
%     xlabel('time');
%     ylabel('Damage factor');
%     hold off
% end
% 
% %% 10. Secondary Observability Analysis
% dt = 0.00001;       
% T_obs = 0.1;    
% T = round(T_obs/dt);
% q = zeros(Nm,N);
% qd = zeros(Nm,N);
% 
% for j = 1:Ndam
%     for i = 1:Nm
%         for k = 1:Nm
%             KKd(i,k) = Kd(j,i,k);
%         end
%     end
%     K_zone{j} = KKd;
% end
% 
% F = Fgf * sin(w_1*(0:N-1)*deltat);
% 
% for n = 2:N
%     qdd = Mgm \ (F(:,n-1) - Kgk*q(:,n-1));
%     qd(:,n) = qd(:,n-1) + deltat*qdd;
%     q(:,n) = q(:,n-1) + deltat*qd(:,n-1);
% end
% 
% for i = 1:Ndam
%     G = eye(Nm)' * K_zone{i} * eye(Nm);      
%     H_at = Phi * (Mgm \ G);         
%     for n = 1:N
%         D_i(:,n) = H_at * q(:,n);       
%     end
%     O(i) = mean(vecnorm(D_i,2,1));      
% end
% 
% for i = 1:2*Nm
%     Hs = [Phi , X ; X , Phi];         
%     for n = 1:N
%         DS_i(:,n) = Hs * [q(:,n); qd(:,n)];     
%     end
%     OS(i) = mean(vecnorm(DS_i,2,1));      
% end
% 
% OS = OS / max(OS);          
% O = O / max(O);          
% 
% %% 11. Final Plotting (Expanded for 12/24 ranges)
% figure(2);
% bar(OS);
% xlabel('No. of modes (disp & Vel: 1-24)');
% ylabel('Normalized observability');
% title('State Observability');
% grid on;
% 
% figure(3);
% bar(O);
% xlabel('Damage zone (1-12)');
% ylabel('Normalized observability');
% title('Zone-wise damage observability');
% grid on;
% 
% % Label Generics for 12/24 sizes
% zone_labels = arrayfun(@(x) sprintf('Z%d', x), 1:12, 'UniformOutput', false);
% mode_labels = arrayfun(@(x) sprintf('m%d', x), 1:24, 'UniformOutput', false);
% 
% figure(4)
% bar(Obsv_fd_s1);
% set(gca, 'XTickLabel', zone_labels);
% title('Observability (Fd) over Time for scheme 1');
% 
% figure(5)
% bar(Obsv_fd_s2);
% set(gca, 'XTickLabel', zone_labels);
% title('Observability (Fd) over Time for scheme 2');
% 
% figure(6)
% bar(Obsv_state_s1);
% set(gca, 'XTickLabel', mode_labels);
% title('Observability Strength (modes) over Time for scheme 1');
% 
% figure(7)
% bar(Obsv_state_s2);
% set(gca, 'XTickLabel', mode_labels);
% title('Observability Strength (modes) over Time for scheme 2');
% 
% figure(8)
% bar(Obsv_ekf(2*Nm+1:2*Nm+Ndam));
% set(gca, 'XTickLabel', zone_labels);
% title('Observability Strength (dams) over Time for EKF');
% 
% P1_s1 = 1./Obsv_state_s1;
% P1_s2 = 1./Obsv_state_s2;
% 
% P2_s1 = 1./Obsv_fd_s1;
% P2_s2 = 1./Obsv_fd_s2;
% %%
% figure(10);
%   for i =1:Nm
%         subplot(2,3,i);
%         plot(linspace(0,deltat*N,N),xx_store_true(i,:),'k-'); hold on;
%         plot(linspace(0,deltat*N,N),xx_store_dkf2(i,:),'m--'); hold on; 
%         plot(linspace(0,deltat*N,N),xx_store_dkf1(i,:),'b--'); hold on;
%         plot(linspace(0,deltat*N,N),xx_store_ekf(i,:),'g--');
%         xlabel('time');
%         ylabel('displacement');
%         title('displacement vs time');
% 
%   end
% % 
% 
% %% 
% % P_ekf = 1./Obsv_ekf;
% % for i=1:24
% %     P_ekf(i)
% % end
% 
% %% cHECKING THE REDUNDANCY 
% % [R, pivots] = rref(Phi);
% % 
% % % 'pivots' will list the columns that are independent.
% % % Any number from 1 to 12 NOT in 'pivots' is a redundant column.
% % redundant_cols = setdiff(1:12, R);
% % disp('Redundant Columns are:');
% % disp(redundant_cols);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%  Lets run the GA 
% %%%%%%%%%%%%%%%5
% %% GA ENABLED EKF FOR BETTER PREDICTION 
% 
% 1. Load your pre-computed data
load('PlateData.mat');

% 2. Define the objective function handle
objective_func = @(x) ekf_objective(x, plate_physics);

% 3. Setup Options with the function handle
try
    options = optimoptions(@ga, ...
        'PopulationSize', 30, ...
        'MaxGenerations', 20, ...
        'Display', 'iter', ...
        'PlotFcn', @gaplotbestf, ...
        'UseParallel', false); % Set to true only if you have Parallel Computing Toolbox

    % 4. Run the GA
    nvars = 5; 
    lb = [-18,-15 ,-18, -18, 0]; 
    ub = [-12,-10, -12,  -12,   2];
    [best_exps, min_val] = ga(objective_func, nvars, [], [], [], [], lb, ub, [], options);

catch ME
    % Friendly error if toolbox is missing
    if strcmp(ME.identifier, 'optimlib:optimoptions:InvalidSolver')
        fprintf('\n--- ERROR: Global Optimization Toolbox NOT FOUND ---\n');
        fprintf('Check your installation by typing: ver\n');
    else
        rethrow(ME);
    end
end

% 5. Display Results
fprintf('\n--- OPTIMIZATION COMPLETE ---\n');
fprintf('Optimized Q_state: 10^(%.2f)\n', best_exps(1));
fprintf('Optimized Q_param: 10^(%.2f)\n', best_exps(2));
fprintf('Optimized R: 10^(%.2f)\n', best_exps(3));


%% GA'S MAIN FITNESS FUNCTION 
function mse_error = ekf_objective(x_exp, plate_physics)
    % x_exp: [logQ_state, logQ_param, logR] - GA's suggested exponents
    N = 5000;
    
    % 1. Reconstruct EKF matrices from exponents
    P_val = 10^x_exp(1);
    Pp_val = 10^x_exp(2);
    Q_val = 10^x_exp(3);
    Qp_val = 10^x_exp(4);
    R_val = 10^x_exp(5);
    
    %% UNPACKING CONSTANTS FROM STRUCTURE 
    Nm = plate_physics.Nm;
    Ndam = plate_physics.Ndam;
    nObserv = plate_physics.nObserv;    
    deltat = plate_physics.deltat;
    w_1 = plate_physics.w_1;
    dim = 2*Nm; % Fixed missing dim definition

    % Numerical Matrices
    Mgm = plate_physics.Mgm;
    Kgun = plate_physics.Kgun;
    Cgc = plate_physics.Cgc;
    Fgf = plate_physics.Fgf;
    Kd = plate_physics.Kd; % The 3D array
    Phi = plate_physics.Phi;
    invMgm = plate_physics.invMgm; 
    
    % Target Data for GA to match
    damf = plate_physics.true_damage; 

    %% INTEGRATION PARAMETERS (Newmark-Beta)
    gamma = 0.5;
    beta = 0.25;
    Haf = Phi; % Assuming Haf maps to your observation matrix Phi
    H = blkdiag(Phi, Phi); % Full state observation matrix

    %% COVARIANCE MATRICES
    % Fixed the data.Nm typo here
    % Restored Pekf matrix so it is defined before the loop
    Pekf = blkdiag(P_val * eye(2*Nm), Pp_val * eye(Ndam));
    Qekf = blkdiag(Q_val * eye(2*Nm), Qp_val * eye(Ndam));
    Rekf = R_val * eye(nObserv);
    
    P = 1e-10 * eye(dim);   
    Pdkf2 = 1e-10 * eye(dim); 
    
    

    Pf = 1e-12*blkdiag(0.3728,0.3122,0.4581,0.3710,0.4272,0.3886);% for 6 sensor scheme 1;
    Pf2 = 1e-12*blkdiag(0.3739 ,0.3122,0.7081, 0.6620,  0.6101 , 0.7058); 
    Q=eye(dim)/10^-14;  %  % scheme 1 -> eye(dim)/10^-10 for 6 sensor (unchange modes)
    Qdkf2=eye(dim)/10^-14; % sche 
    Qf = Pf; 
    Qf2 = Pf2; 

    R=1*eye(2*nObserv)*10^(4);  % for state ( disp & vel) of Scheme - 1 for 6 sen ->1*eye(2*nObserv)*10^(-4) works better ( no noise)  
    Rf=1*eye(nObserv)*10^(2);
    Rf2=1*eye(nObserv)*10^(2);
    Rdkf2=1*eye(2*nObserv)*10^(4);
    
    %% INITIALIZE ALL STATE VARIABLES (Crucial to prevent crashes)
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

gamma = 1/2;
beta = 1/2;
x_st = zeros(Nm,1);
v_st = zeros(Nm,1);
acc = zeros(Nm,1);
ac = zeros(Nm,1);
Fg = zeros(Nm,1);


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

    v_st(:) = vp + deltat*(1-gamma)*accpri + deltat*gamma*accri;
    x_st(:) = xp + deltat*vp + (1-2*beta)/2*deltat^2*accpri + beta*deltat^2*accri;
    accri(1) = ac(1) + 1e-4*(rand-0.5)/5;
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
    vv = (xx(1:Nm) - xxp(1:Nm)) / deltat;
    aa = (vv(1:Nm) - vvp(1:Nm)) / deltat;
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


    %% 3. Calculate Fitness Score
    % Compare EKF prediction directly to True Damage
    diff = damf(:, 1:N) - fd_est_ex(:, 1:N);
    mse_error = mean(diff(:).^2); 
    
    % Penalty for divergence
    if isnan(mse_error) || mse_error > 1e10
        mse_error = 1e10; 
    end
end
