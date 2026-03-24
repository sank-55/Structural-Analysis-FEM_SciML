% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %    %% 6 Modes, 12 Damage Factors 10 SENSOR AT MOST OBS LOCATION  - Dual KF & EKF with Observability
% %      ONLY EKF HAVE TO SEE ( REST NOISE , P,W IS ADJUSTED ) 
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
nObserv = 10;     %(8 to 12 easily predicted by 6 MODES )     % No of observable zones (sensors)

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
Pekf = 1e4*blkdiag(1.3421e-14,1.2625e-16, 7.7437e-16, 4.8198e-17,1.0609e-17,7.7538e-17,5.5859e-09, 4.6491e-11,2.8968e-10 ,1.7384e-11, 5.2969e-12, 2.8044e-11, 1.9752e-08,6.7107e-09,6.8089e-09,2.0973e-08,2.2424e-08,3.8958e-09,3.9001e-09, 2.3562e-08,1.4734e-08,9.9054e-09 ,9.9128e-09,1.4875e-08 );

Qekf = 1e-8 * eye(2*Nm + Ndam);
Qf = Pf; 
Qf2 = Pf2; 

R = 10^-1 * eye(2*nObserv);  
Rf = 10^-0.2 * eye(nObserv);
Rdkf2 = 10^-1 * eye(2*nObserv); 
Rf2 = 10^-0.2 * eye(nObserv);
Rekf = 10^-0 * eye(nObserv);

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
