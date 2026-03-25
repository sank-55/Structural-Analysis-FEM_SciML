%% MAIN SCRIPT: GA Optimization for EKF Damage Detection
clear; clc;

fprintf('--- STEP 1: Pre-computing Symbolic Matrices (Runs ONCE) ---\n');
fprintf('Please wait, integrating mode shapes... ');

% 1. Plate Parameters
syms x y
p = 1000; 
ro = 2700;           % Density (kg/m3)
yo = 70*10^9;        % Young's modulus
neu = 0.3;           % Poisson's ratio
h = 0.001;           % Thickness (m)
a = 0.6;             % Length (m)
b = 0.4;             % Breadth (m)
D = (yo * h^3) / (12 * (1 - neu^2));

Nm = 6;              % 6 Modes used in definition
Nn = Nm;
Ndam = 12;           % 12 Damage zones
N = 5000;
deltat = 0.0001;
w_1 = 100*pi;

% 2. Mode Shapes & Derivatives
phi = [sin(1*pi*x/a)*sin(1*pi*y/b), sin(1*pi*x/a)*sin(2*pi*y/b), ...
       sin(2*pi*x/a)*sin(1*pi*y/b), sin(2*pi*x/a)*sin(2*pi*y/b), ...
       sin(1*pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(1*pi*y/b)];

dxxphi = diff(phi, x, 2);
dxyphi = diff(diff(phi, x), y);
dyyphi = diff(phi, y, 2);

% 3. Coordinates
damdox = [0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a; 0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a; 0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a];
damdoy = [0 b/3; 0 b/3; 0 b/3; 0 b/3; b/3 2*b/3; b/3 2*b/3; b/3 2*b/3; b/3 2*b/3; 2*b/3 b; 2*b/3 b; 2*b/3 b; 2*b/3 b];

sensox = [ 0.1220  0.4780  0.1220  0.4780  0.4881  0.1119  0.2949  0.3051  0.3966  0.2034   ];
sensoy = [ 0.0718  0.0718  0.3282  0.3282  0.2359  0.2359  0.3282  0.0718  0.1538  0.1538   ];
nObserv = 10; % Dynamically set to 12

% 4. Matrix Pre-allocation & Symbolic Integration
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

for i = 1:Ndam
    for j = 1:Nm
        for k = 1:Nm
            Kd(i,j,k) = double(D * int(int(dxxphi(j)*dxxphi(k) + 2*dxyphi(j)*dxyphi(k) + dyyphi(j)*dyyphi(k), x, damdox(i,1), damdox(i,2)), y, damdoy(i,1), damdoy(i,2)));
        end
    end
end

Cgc = 0.0003 * Mgm + 0.0003 * Kgk;     
Kgun = Kgk;  

% 5. Observation Matrices
Phi = zeros(nObserv, Nm);
for i = 1:nObserv
    x1 = sensox(i); y1 = sensoy(i);
    Phi(i,:) = [sin(1*pi*x1/a)*sin(1*pi*y1/b), sin(1*pi*x1/a)*sin(2*pi*y1/b), ...
                sin(2*pi*x1/a)*sin(1*pi*y1/b), sin(2*pi*x1/a)*sin(2*pi*y1/b), ...
                sin(1*pi*x1/a)*sin(3*pi*y1/b), sin(3*pi*x1/a)*sin(1*pi*y1/b)];
end
X_mat = zeros(nObserv, Nm);
H = [Phi X_mat ; X_mat Phi];
Haf = Phi;

% 6. True Damage Scenario
damf = zeros(Ndam, N);
for i = 200:N,  damf(1,i) = (i-200)^.7 * 0.3 / (N-200)^.7; end
for i = 1200:N, damf(2,i) = (i-1200) * 0.6 / (N-1200); end
for i = 600:N,  damf(3,i) = (i-600) * 0.3 / (N-600); end
for i = 30:N,   damf(5,i) = (i-30) * 0.2 / (N-30); end
for i = 800:N,  damf(8,i) = (i-800) * 0.4 / (N-800); end
for i = 1500:N, damf(11,i) = (i-1500) * 0.25 / (N-1500); end

% 7. Pack into Struct for GA
plate.Mgm = Mgm; plate.Kgk = Kgk; plate.Fgf = Fgf; plate.Kd = Kd; 
plate.Cgc = Cgc; plate.Kgun = Kgun; plate.Phi = Phi; plate.H = H; 
plate.Haf = Haf; plate.damf = damf; plate.Nm = Nm; plate.Ndam = Ndam; 
plate.N = N; plate.deltat = deltat; plate.w_1 = w_1; plate.nObserv = nObserv;
fprintf('Done!\n');

%% --- STEP 2: Run Genetic Algorithm ---
objective_func = @(x) ekf_objective(x, plate);

try
    fprintf('Starting GA Optimization...\n');
    options = optimoptions(@ga, 'PopulationSize', 50, 'MaxGenerations', 50, ...
        'Display', 'iter', 'PlotFcn', @gaplotbestf, 'UseParallel', false); 

    nvars = 5; 
    lb = [-16,-16, -12, -12, 1]; 
    ub = [-8, -8,  -8,  -1,   4];
    [best_exps, min_val] = ga(objective_func, nvars, [], [], [], [], lb, ub, [], options);

    fprintf('\n--- OPTIMIZATION COMPLETE ---\n');
    
    fprintf('Optimized P_state: 10^(%.2f)\n', best_exps(1));
    fprintf('Optimized p_param: 10^(%.2f)\n', best_exps(2));
    fprintf('Optimized Q_state: 10^(%.2f)\n', best_exps(3));
    fprintf('Optimized Q_param: 10^(%.2f)\n', best_exps(4));
    fprintf('Optimized R: 10^(%.2f)\n', best_exps(3));
catch ME
    rethrow(ME);
end

%% --- STEP 3: The Objective Function (Cleaned & Optimized) ---
function mse_error = ekf_objective(x_exp, plate)
    % Unpack parameters
    P_val = 10^x_exp(1); Pp_val = 10^x_exp(2); R_val = 10^x_exp(5);
    Q_val = 10^x_exp(3);Qp_val = 10^x_exp(4);
    % Unpack structural matrices
    Mgm = plate.Mgm; Kgk = plate.Kgk; Fgf = plate.Fgf; Kd = plate.Kd; 
    Cgc = plate.Cgc; Kgun = plate.Kgun; Phi = plate.Phi; H = plate.H; 
    Haf = plate.Haf; damf = plate.damf; Nm = plate.Nm; Ndam = plate.Ndam; 
    N = plate.N; deltat = plate.deltat; w_1 = plate.w_1; nObserv = plate.nObserv;
    dim = 2 * Nm;

    %% Filter Tuning Matrices
    P = 1e-10 * eye(dim);   
    Pdkf2 = 1e-10 * eye(dim); 
    Pf = 1e-8 * blkdiag(0.1246, 0.0655,0.0535, 0.1486,0.1804,0.0217, 0.0227,0.1692,0.1699,0.081, 0.0808, 0.1641); 
    Pf2 = 1e-8 * blkdiag(0.1246, 0.0655,0.0535, 0.1486,0.1804,0.0217, 0.0227,0.1692,0.1699,0.081, 0.0808, 0.1641);

    Q = eye(dim) / 1e-10;  
    Qdkf2 = eye(dim) / 1e-10; 
    %Pekf = 1e4*blkdiag(1.3421e-14,1.2625e-16, 7.7437e-16, 4.8198e-17,1.0609e-17,7.7538e-17,5.5859e-09, 4.6491e-11,2.8968e-10 ,1.7384e-11, 5.2969e-12, 2.8044e-11, 1.9752e-08,6.7107e-09,6.8089e-09,2.0973e-08,2.2424e-08,3.8958e-09,3.9001e-09, 2.3562e-08,1.4734e-08,9.9054e-09 ,9.9128e-09,1.4875e-08 );
   % Pekf = 1e5*blkdiag(1.3421e-14,1.2625e-16, 7.7437e-16, 4.8198e-17,1.0609e-17,7.7538e-17,5.5859e-09, 4.6491e-11,2.8968e-10 ,1.7384e-11, 5.2969e-12, 2.8044e-11, 1.9752e-08,6.7107e-09,6.8089e-09,2.0973e-08,2.2424e-08,3.8958e-09,3.9001e-09, 2.3562e-08,1.4734e-08,9.9054e-09 ,9.9128e-09,1.4875e-08 );
    Pekf = blkdiag(P_val * eye(2*Nm), Pp_val * eye(Ndam)); 
    Qekf = blkdiag(Q_val * eye(2*Nm), Qp_val * eye(Ndam));
    Rekf = R_val * eye(nObserv);
    Qf = Pf; Qf2 = Pf2; 

    R = 10^-1 * eye(2*nObserv);  
    Rf = 10^-0.2 * eye(nObserv);
    Rdkf2 = 10^-1 * eye(2*nObserv); 
    Rf2 = 10^-0.2 * eye(nObserv);

    %% Variable Initialization
    xx = zeros(dim, 1);
    acc = zeros(Nm, 1); ac = zeros(Nm, 1);
    fac = zeros(Ndam, 1); facc = zeros(Ndam, 1);
    fd_dkf2 = zeros(Ndam, 1);
    fd_est_ex = zeros(Ndam, N);
    
    xxp = zeros(dim, 1); xxp2 = zeros(dim, 1); xxp_exfd = zeros(dim+Ndam, 1);  
    vp = zeros(Nm, 1); xp = zeros(Nm, 1);
    accp = zeros(Nm, 1); acp = zeros(Nm, 1); acp_2 = zeros(Nm, 1);
    xrip = zeros(Nm, 1); vrip = zeros(Nm, 1); accpri = zeros(Nm, 1);
    x_st = zeros(Nm,1); v_st = zeros(Nm,1);
    
    xi = zeros(Nm, Ndam);
    gamma_val = 0.5; beta_val = 0.5;

    %% Main Time Integration Loop
    for i = 2:N
        % True System Update
        Kgk_true = Kgun;
        for j = 1:Ndam
            tem = squeeze(Kd(j,:,:)); % Vectorized speedup
            Kgk_true = Kgk_true - damf(j,i) * tem;
        end

        AA = Mgm + Cgc*deltat*gamma_val + Kgk_true*beta_val*deltat^2;
        BB = Cgc*deltat*(1-gamma_val) + Kgk_true*(1-2*beta_val)*deltat^2/2;
        CC = Cgc + deltat*Kgk_true;

        t = (i-1)*deltat;
        Fg = Fgf * sin(w_1*t);
        acc(:) = AA \ (Fg - BB*accp - CC*v_st - Kgk_true*x_st);
        accp = acc;

        for kk = 1:Nm
            ac(kk) = acc(kk) + 1e-4*(randn-0.5)/50;
        end
        accri = ac + 1e-4*(randn-0.5)/70;
        accri(1) = ac(1) + 1e-4*(rand-0.5)/5;

        v_st(:) = vp + deltat*(1-gamma_val)*accpri + deltat*gamma_val*accri;
        x_st(:) = xp + deltat*vp + (1-2*beta_val)/2*deltat^2*accpri + beta_val*deltat^2*accri;
        vp = v_st; xp = x_st;

        vri = vrip + (acp+ac)/2 * deltat;
        xri = xrip + (vrip+vri)/2 * deltat;

        % FIXED DIMENSIONS: Observation Vectors
        z_pos = Phi * (xri + xri .* (0.02 * rand(Nm,1)));
        z_vel = Phi * (vri + vri .* (0.02 * rand(Nm,1)));
        z = [z_pos; z_vel]; 
        zf = Phi * (ac + ac .* (0.03 * rand(Nm,1)));

        %% Filter Prediction
        disp_est = xxp(1:Nm) + xxp(Nm+1:2*Nm)*deltat;
        vel_est = xxp(Nm+1:2*Nm) + acp(1:Nm)*deltat;
        fac = facc;
        Kdx = zeros(Nm, Ndam);

        for j = 1:Ndam
            tem = squeeze(Kd(j,:,:));
            xi(:,j) = tem * disp_est;
            Kdx(:,j) = xi(:,j);
        end

        F = Fgf * sin(w_1*t);  
        Af = eye(Ndam);
        Hf = Haf * (Mgm \ Kdx);
        Cgcx = Mgm \ (Kgun*disp_est + Cgc*vel_est);
        HBf = Haf * ((Mgm \ F) - Cgcx);

        % Scheme 2 Stiffness
        Kgks2 = Kgun;
        for j = 1:Ndam
            tem = squeeze(Kd(j,:,:));
            Kgks2 = Kgks2 - fd_dkf2(j) * tem; 
        end

        Kdx2 = zeros(Nm, Ndam);
        for j = 1:Ndam
            tem = squeeze(Kd(j,:,:));
            xi(:,j) = tem * (xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2);
            Kdx2(:,j) = xi(:,j);
        end

        InvMatS2 = Mgm + deltat*0.5*Cgc + deltat^2*0.25*Kgks2;
        Hdkf2 = Phi * (InvMatS2 \ Kdx2);
        a_dkf2 = Hdkf2 * fd_dkf2 + Phi * (InvMatS2 \ (F - Kgun*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*0.25*acp_2) - Cgc*(xxp2(Nm+1:2*Nm) + deltat^2*0.5*acp_2)));

        % Update Damage Factors
        Pnf = Af * Pf * Af' + Qf;
        Kf = Pnf * Hf' / ((Hf * Pnf * Hf') + Rf);
        fac = facc + Kf * (zf - Hf*fac(:) - HBf); % Fixed dimension

        Pnf2 = Af * Pf2 * Af' + Qf2;
        Kdkf2 = Pnf2 * Hdkf2' / ((Hdkf2 * Pnf2 * Hdkf2') + Rf2);
        fd_dkf2 = fd_dkf2 + Kdkf2 * (zf - a_dkf2); % Fixed dimension

        %% 2nd Loop for DEKF
        Kgk1 = Kgun;
        for j = 1:Ndam
            tem = squeeze(Kd(j,:,:));
            Kgk1 = Kgk1 - fac(j)*tem;  
        end

        Ac = [zeros(Nm,Nm) eye(Nm,Nm); -(Mgm \ Kgk1) -(Mgm \ Cgc)];
        Bc = [zeros(Nm,Nm); Mgm \ eye(Nm)];

        Kgdkf2 = Kgun;
        for j = 1:Ndam
            tem = squeeze(Kd(j,:,:));
            Kgdkf2 = Kgdkf2 - fd_dkf2(j)*tem;  
        end

        Adkf2 = [zeros(Nm,Nm) eye(Nm,Nm); -(Mgm \ Kgdkf2) -(Mgm \ Cgc)];
        Bdkf2 = [zeros(Nm,Nm); Mgm \ eye(Nm)];

        S_ekf = zeros(Nm, Ndam);
        for j = 1:Ndam
           tem = squeeze(Kd(j,:,:));
           S_ekf(:,j) = tem * xxp_exfd(1:Nm);
        end

        Aex = [zeros(Nm,Nm) , eye(Nm,Nm) , zeros(Nm,Ndam);
               -(Mgm \ Kgk1) , -(Mgm \ Cgc) , (Mgm \ S_ekf);
               zeros(Ndam,Nm), zeros(Ndam,Nm), eye(Ndam,Ndam)];
        Bex = [zeros(Nm,Nm); Mgm \ eye(Nm); zeros(Ndam,Nm)];

        % Discretization (Euler Expansion)
        A = eye(dim); BB = eye(dim);
        Ax2 = eye(dim); Bx2_ = eye(dim);
        Ax = eye(dim+Ndam); BBx = eye(dim+Ndam);
        
        for jj = 1:10
            A = A + (Ac*deltat)^jj / factorial(jj);
            BB = BB + (Ac*deltat)^jj / factorial(jj+1);
            Ax2 = Ax2 + (Adkf2*deltat)^jj / factorial(jj);
            Bx2_ = Bx2_ + (Adkf2*deltat)^jj / factorial(jj+1);
            Ax = Ax + (Aex*deltat)^jj / factorial(jj);
            BBx = BBx + (Aex*deltat)^jj / factorial(jj+1);
        end
        B = BB * Bc * deltat;
        Bx2 = Bx2_ * Bdkf2 * deltat;
        Bx = BBx * Bex * deltat;

        % EKF Update
        xx_exfd = Ax * xxp_exfd + Bx * F;
        Pekfn = Ax * Pekf * Ax' + Qekf; 

        Kgk2 = Kgun;
        for j = 1:Ndam
            tem = squeeze(Kd(j,:,:));
            Kgk2 = Kgk2 - xx_exfd(2*Nm+j) * tem;  
        end

        for j = 1:Ndam
           tem = squeeze(Kd(j,:,:));
           S_ekf(:,j) = tem * xx_exfd(1:Nm);
        end

        Hekf = [-Phi*(Mgm \ Kgk2) , -Phi*(Mgm \ Cgc) , Phi*(Mgm \ S_ekf) ];
        z_diff_ekf = zf - Phi*(Mgm \ (F - Kgk2*xx_exfd(1:Nm) - Cgc*xx_exfd(Nm+1:2*Nm)));
        Kekf = Pekfn * Hekf' / ((Hekf * Pekfn * Hekf') + Rekf);
        xxp_exfd = xx_exfd + Kekf * z_diff_ekf;

        fd_est_ex(:,i) = xxp_exfd(2*Nm +1: end);

        % Scheme updates 
        xx = A*xxp + B*F;
        Pn = A*P*A' + Q;    
        K = Pn * H' / ((H*Pn*H') + R);
        xx = xx + K*(z - H*xx); % Fixed Dimension   

        xx2 = Ax2*xxp2 + Bx2*F;
        Pndkf2 = Ax2*Pdkf2*Ax2' + Qdkf2;    
        Kx2 = Pndkf2 * H' / ((H*Pndkf2*H') + Rdkf2);
        xx2 = xx2 + Kx2*(z - H*xx2); % Fixed Dimension  

        % Covariances & Stores
        facc = fac;
        P = (eye(dim) - K*H) * Pn;
        Pf = (eye(Ndam) - Kf*Hf) * Pnf;
        Pdkf2 = (eye(dim) - Kx2*H) * Pndkf2;
        Pf2 = (eye(Ndam) - Kdkf2*Hdkf2) * Pnf2;

        xxp = xx; acp_2 = (xx2(Nm+1:2*Nm) - xxp2(Nm+1:2*Nm)) / deltat;
        xxp2 = xx2; accp = acc; acp = ac;
        xrip = xri; vrip = vri; accpri = accri;
    end
    
    %% 3. Calculate Fitness Score
    error_diff = damf(:, 1:N) - fd_est_ex(:, 1:N);
    mse_error = mean(error_diff(:).^2); 
    
    if isnan(mse_error) || mse_error > 1e10 || ~isreal(mse_error)
        mse_error = 1e10; 
    end
end
