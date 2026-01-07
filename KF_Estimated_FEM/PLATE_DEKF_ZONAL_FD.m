%% ==========================================================
%   RECOND , Jump in the estimated data Concept is working , Have to look into this

%% ------------------------------------------------------------------
clear; close all; clc;

%% ---------------- geometry & material ----------------
E = 70e9;         % Young's modulus (Pa)
h = 0.005;         % Plate thickness (m)
rho = 7050;        % Density (kg/m³)
nu = 0.3;          % Poisson's ratio

Lx = 1.2; Ly = 0.4;    % Plate dimensions (m)
q0 = 1.5;        % Distributed load (N/m²)

% Mesh configuration
nx = 24; ny = 8;  % balanced mesh
dx = Lx/nx; dy = Ly/ny;
Nelem = nx*ny;
Nnode = (nx+1)*(ny+1);
Ndof = 3*Nnode;    % [w, theta_x, theta_y] per node

fprintf('=== PLATE FEM WITH DKF (fixed) ===\n');
fprintf('Mesh: %dx%d, Nodes: %d, Elements: %d\n', nx, ny, Nnode, Nelem);

%% ---------------- BENDING D MATRIX ----------------
D = (E*h^3 / (12*(1 - nu^2))) * [1, nu, 0;
                                 nu, 1, 0;
                                 0, 0, (1 - nu)/2];

%% ---------------- NODES & CONNECTIVITY ----------------
coords = zeros(Nnode,2);
ncount = 1;
for j = 0:ny
    for i = 0:nx
        coords(ncount,:) = [i*dx, j*dy];
        ncount = ncount + 1;
    end
end

% Element connectivity (quad)
c = zeros(Nelem,4);
elem = 0;
for j = 1:ny
    for i = 1:nx
        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n3 = n2 + (nx+1);
        n4 = n1 + (nx+1);
        elem = elem + 1;
        c(elem,:) = [n1 n2 n3 n4];
    end
end

%% ---- take the centroid and boundary of the zones --------
elemCenter = zeros(Nelem,2 ); % store element centroid

%% ---------------- SHAPE FUNCTIONS & GAUSS INTEGRATION ----------------
gp = [-1/sqrt(3), 1/sqrt(3)];
wgp = [1,1];

syms xi eta real
% Hermite cubic shape functions (as in your original)
H1 = 1 - 3*xi^2 + 2*xi^3;
H2 = (xi - 2*xi^2 + xi^3);
H3 = 3*xi^2 - 2*xi^3;
H4 = (-xi^2 + xi^3);
H_xi = [H1; H2; H3; H4];
H_eta = subs(H_xi, xi, eta);

% 12 shape functions for Kirchhoff plate
N1  = H_xi(1) * H_eta(1); N2  = H_xi(2) * H_eta(1);
N3  = H_xi(1) * H_eta(2); N4  = H_xi(3) * H_eta(1);
N5  = H_xi(4) * H_eta(1); N6  = H_xi(3) * H_eta(2);
N7  = H_xi(3) * H_eta(3); N8  = H_xi(4) * H_eta(3);
N9  = H_xi(3) * H_eta(4); N10 = H_xi(1) * H_eta(3);
N11 = H_xi(2) * H_eta(3); N12 = H_xi(1) * H_eta(4);

N_sym = [N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12];
N_xixi  = simplify(diff(N_sym, xi, 2));
N_etaeta = simplify(diff(N_sym, eta, 2));
N_xieta = simplify(diff(diff(N_sym, xi), eta));

detJ_const = (2/dx)*(2/dy);

% Precompute at Gauss points
numGP = 4;
N_at_gp = zeros(numGP, 12);
B_at_gp = zeros(3,12,numGP);
wts_at_gp = zeros(numGP,1);
cnt = 0;
for i = 1:2
    for j = 1:2
        cnt = cnt + 1;
        xi_val = gp(i); eta_val = gp(j);
        Nval = double(subs(N_sym, [xi,eta], [xi_val, eta_val]));
        N_xixi_val = double(subs(N_xixi, [xi,eta], [xi_val, eta_val]));
        N_etaeta_val = double(subs(N_etaeta, [xi,eta], [xi_val, eta_val]));
        N_xieta_val = double(subs(N_xieta, [xi,eta], [xi_val, eta_val]));
        Bval = (4/dx^2) * [N_xixi_val; N_etaeta_val; 2*N_xieta_val];
        N_at_gp(cnt,:) = Nval;
        B_at_gp(:,:,cnt) = Bval;
        wts_at_gp(cnt) = wgp(i)*wgp(j);
    end
end

%% ---------------- ASSEMBLY ----------------
Ke_store = cell(Nelem,1); % store element stiffness
elemDof = cell(Nelem,1);

fprintf('Assembling global matrices...\n');
Kg = sparse(Ndof, Ndof);
Mg = sparse(Ndof, Ndof);
Fg = zeros(Ndof,1);

for e = 1:Nelem
    nodes = c(e,:);
    % store centroid
    elemCenter(e,:) = mean(coords(nodes,:),1);

    Ke = zeros(12,12);
    Me = zeros(12,12);
    Fe = zeros(12,1);
    for gpidx = 1:numGP
        Nval = N_at_gp(gpidx,:)';
        Bval = B_at_gp(:,:,gpidx);
        wt = wts_at_gp(gpidx);
        Ke = Ke + (Bval')*D*Bval*detJ_const*wt;
        Me = Me + (rho*h) * (Nval * Nval') * detJ_const * wt;
        Fe = Fe + q0 * detJ_const * wt * Nval;
    end
    dofs = reshape([3*nodes-2; 3*nodes-1; 3*nodes],1,[]);
    elemDof{e} = dofs;
    Ke_store{e} = Ke;
    Kg(dofs,dofs) = Kg(dofs,dofs) + Ke;
    Mg(dofs,dofs) = Mg(dofs,dofs) + Me;
    Fg(dofs) = Fg(dofs) + Fe;
end

%% ------------Define zones (geometric split in x) ------
Nz = 3;
zoneID = zeros(Nelem,1);
for e = 1:Nelem
    x = elemCenter(e,1);
    if x <= Lx/3 %&& y <= Ly/2
        zoneID(e) = 1;
    elseif x <= 2*Lx/3
        zoneID(e) = 2;
    else
        zoneID(e) = 3;
    end
end

% zone centers
zone_center = zeros(Nz,2);
for z = 1:Nz
    elems = find(zoneID==z);
    zone_center(z,:) = mean(elemCenter(elems,:),1);
end

% Build zone stiffness matrices (full DOF)
Kz = cell(Nz,1);
for z = 1:Nz
    Kz{z} = sparse(Ndof, Ndof);
end
for e = 1:Nelem
    z = zoneID(e);
    dofs = elemDof{e};
    Ke = Ke_store{e};
    Kz{z}(dofs,dofs) = Kz{z}(dofs,dofs) + Ke;
end

% check zoning assembly
Ksum = sparse(Ndof,Ndof);
for z=1:Nz, Ksum = Ksum + Kz{z}; end
err = norm(Ksum - Kg, 'fro') / norm(Kg, 'fro');
fprintf('Zoning error = %.2e\n', err);

%% ---------------- BOUNDARY CONDITIONS ----------------
tol = 1e-9;
left_nodes = find(abs(coords(:,1)) < tol);   % left edge clamped
fixedDOF = reshape([3*left_nodes-2; 3*left_nodes-1; 3*left_nodes],1,[]);
freeDOF = setdiff(1:Ndof, fixedDOF);

Kf = Kg(freeDOF, freeDOF);
Kf_undam = Kf;
Mf = Mg(freeDOF, freeDOF);
Ff = Fg(freeDOF);

% apply BC to zone stiffness
Kzf = cell(Nz,1);
for z = 1:Nz
    Kzf{z} = Kz{z}(freeDOF, freeDOF);
end

fprintf('Free DOF: %d, Fixed DOF: %d\n', length(freeDOF), length(fixedDOF));

%% ----- modal analysis ---------------
nModesKeep = 6;     % keep 10 modes (adjust if necessary)
opts.isreal = 1;
opts.issym = 1;
[PHI_all, W2_all] = eigs(Kf_undam, Mf, nModesKeep, 'smallestabs', opts);
omega_all = sqrt(abs(diag(W2_all)));

% Construct Rayleigh damping using two modes
zeta_target = 0.02;
w1 = omega_all(1); w2 = omega_all(min(3,length(omega_all)));
A = [1/w1, w1; 1/w2, w2]; b = 2*zeta_target*[1;1];
alpha_beta = A\b; alphaR = alpha_beta(1); betaR = alpha_beta(2);
Cf = alphaR * Mf + betaR * Kf_undam;

modes = nModesKeep;
Phi_modal = PHI_all(:,1:modes);

% modal reduced matrices
K_modal_red = Phi_modal' * Kf_undam * Phi_modal;
M_modal_red = Phi_modal' * Mf * Phi_modal;
C_modal_red = Phi_modal' * Cf * Phi_modal;
F_modal_red = Phi_modal' * Ff;

% zone stiffness in modal space
Kz_modal_red = cell(Nz,1);
for z = 1:Nz
    Kz_modal_red{z} = Phi_modal' * Kzf{z} * Phi_modal;
end

%% -------- assigning sensors → center of each zone (nearest node) ----
sensorNode = zeros(Nz,1);
for z = 1:Nz
    dist = vecnorm(coords - zone_center(z,:), 2, 2);
    [~, sensorNode(z)] = min(dist);
end

% show mesh + zones + sensors
figure('Name','Zones & Sensors'); hold on; axis equal;
cmap = lines(Nz);
for e = 1:Nelem
    nodes = c(e,:);
    X = coords(nodes,1);
    Y = coords(nodes,2);
    z = zoneID(e);
    patch(X, Y, cmap(z,:), 'EdgeColor','k');
end
for z = 1:Nz
    elems = find(zoneID==z);
    nodes = unique(c(elems,:));
    xc = mean(coords(nodes,1)); yc = mean(coords(nodes,2));
    text(xc, yc, num2str(z),'FontWeight','bold');
end
plot(coords(sensorNode,1), coords(sensorNode,2), 'rp','MarkerSize',12,'MarkerFaceColor','r');
title('Zones and Sensors'); xlabel('x'); ylabel('y');
drawnow;

%% ---- Set realistic time stepping for true response (avoid tiny dt) ----
dt = 1e-53;              % time step (s) — choose based on modal highest frequency
Ttotal = 1e-50;           % total time in seconds
tvec = 0:dt:Ttotal;
nSteps = length(tvec)-1;
fprintf('Time steps: %d, dt = %.2e s\n', nSteps, dt);

%% ---- define true damage time-history for synthetic truth ----
fd = zeros(Nz, nSteps+1);
% Example: zone 3 ramps to 0.4, zone2 ramps to 0.2, zone1 healthy
for z = 1:Nz
    if z==3
        fd(z,:) = linspace(0,0.4,nSteps+1);
    elseif z==2
        fd(z,:) = linspace(0,0.2,nSteps+1);
    else
        fd(z,:) = zeros(1,nSteps+1);
    end
end

%% -- initialize for Newmark true solution (full free DOFs)
nFree = size(Kf,1);
u_cur = zeros(nFree,1);   % displacement at free DOFs
v_cur = zeros(nFree,1);   % velocity
a_cur = zeros(nFree,1);   % acceleration

% preallocate stores (transverse w only)
w_dofs = 3*(1:Nnode) - 2;            % global w DOF indices
w_dofs_free = intersect(w_dofs, freeDOF);
Nwfree = length(w_dofs_free);

u_store_true_w = zeros(Nnode, nSteps+1);
v_store_true_w = zeros(Nnode, nSteps+1);
a_store_true_w = zeros(Nnode, nSteps+1);

% Newmark constants
beta_nm = 0.25; gamma_nm = 0.5;

% Pre-factor (we assemble Kf_dam each step)
fprintf('Simulating true response (Newmark)...\n');
for step = 1:nSteps
    t = tvec(step);
    % assemble damaged stiffness for this step (reset each time)
    Kf_dam = sparse(nFree, nFree);
    for z = 1:Nz
        Kf_dam = Kf_dam + (1 - fd(z,step)) * Kzf{z};
    end

    % effective mass-stiff system for Newmark
    a_cur = Mf \ (Ff - Kf_dam*u_cur - Cf*v_cur);   % initial accel
    u_pred = u_cur + dt*v_cur + (0.5 - beta_nm)*dt^2 * a_cur;
    v_pred = v_cur + (1 - gamma_nm)*dt * a_cur;

    % recompute stiffness at next substep (use fd at step+1)
    Kf_dam_next = sparse(nFree, nFree);
    for z = 1:Nz
        Kf_dam_next = Kf_dam_next + (1 - fd(z,step+1)) * Kzf{z};
    end

    a_next = Mf \ (Ff - Kf_dam_next*u_pred - Cf*v_pred);
    v_next = v_pred + gamma_nm*dt*a_next;
    u_next = u_pred + beta_nm*dt^2 * a_next;

    % store transverse DOFs (rebuild full U)
    Ufull = zeros(Ndof,1);
    Ufull(freeDOF) = u_next;
    w_s = Ufull(w_dofs);
    u_store_true_w(:, step+1) = w_s;

    Ufull = zeros(Ndof,1);
    Ufull(freeDOF) = v_next;
    v_store_true_w(:, step+1) = Ufull(w_dofs);

    Ufull = zeros(Ndof,1);
    Ufull(freeDOF) = a_next;
    a_store_true_w(:, step+1) = Ufull(w_dofs);

    % shift
    u_cur = u_next; v_cur = v_next; a_cur = a_next;
    if mod(step,round(nSteps/10))==0
        fprintf('  step %d/%d t=%.4f s\n', step, nSteps, t);
    end
end
% ensure first column (initial) contain initial values (already zeros)

%% ---- Build mapping Phi_w and sensor selection S ----
% rows in Phi_modal correspond to free DOFs only (size nFree x modes)
% we need rows corresponding to free w DOFs
w_is = ismember(freeDOF, w_dofs);        % logical index into freeDOF for all w DOFs
Phi_w = Phi_modal(w_is, :);              % Nwfree x modes

% sensor DOFs: global w DOFs of sensor nodes
sensor_w_global = 3*sensorNode - 2;      % size Nz x 1 (global DOF index)
% convert to index among free w DOFs
[found, sensor_idx] = ismember(sensor_w_global, w_dofs_free);
if ~all(found)
    error('One or more sensor nodes are constrained (on boundary). Move sensor.');
end

Ns = length(sensor_idx);
S = zeros(Ns, Nwfree);
for i = 1:Ns
    S(i, sensor_idx(i)) = 1;
end

% sensor-mode matrix (maps modal coords to sensor w)
Phi_sen = S * Phi_w;   % size Ns x modes

%% ------- DKF: set covariances & storage ----------------
params = Nz;

% reasonable noise levels (tune for your problem)
a_noise_std = 1e-4;      % acceleration sensor noise (m/s^2)
u_noise_std = 1e-4;      % small displacement noise for u_meas (if used)
v_noise_std = 1e-4;

% Covariances
P_x = 1e-4 * eye(2*modes);         % state covariance
P_fd = 1e-3 * eye(params);         % parameter covariance

Q_x = 1e-5 * eye(2*modes);
Q_fd = 1e-4 * eye(params);

R_u = 1.7 * eye(Ns);
R_v = 1.7 * eye(Ns);
R_x = blkdiag(R_u, R_v);           % measurement covariance for u and v (2*Ns x 2*Ns)
R_fd = 1.7 * eye(Ns);  % accel measurement covariance for parameters

% initial estimates
fd_cur = zeros(params,1);          % start with healthy guess
x_cur = zeros(2*modes,1);          % modal states (q; qdot)

% storage
x_state = zeros(2*modes, nSteps+1);
fd_est = zeros(params, nSteps+1);

u_store_dkf_w = zeros(Nnode, nSteps+1);
a_store_dkf_w = zeros(Nnode, nSteps+1);

% initial measured acceleration at sensors (from true simulation) + noise
a_meas = zeros(Ns, nSteps+1);
for i = 1:Ns
    a_meas(i,1) = a_store_true_w(sensorNode(i),1) + a_noise_std*randn;
end

fprintf('Running DKF...\n');
% Pre-allocate for speed
Istates = eye(2*modes);
Iparams = eye(params);

% main DKF loop
for k = 1:nSteps
    %% --- measurement from "sensors" (synthetic, with noise) ---
    % true acceleration at sensors (use true simulation) + noise
    for i = 1:Ns
        a_meas(i,k+1) = a_store_true_w(sensorNode(i), k+1) + a_noise_std*randn;
    end
    % also displacement and velocity measurements at sensors if needed
    u_meas = zeros(Ns,1);
    v_meas = zeros(Ns,1);
    for i = 1:Ns
        u_meas(i) = u_store_true_w(sensorNode(i), k+1) + u_noise_std*randn;
        v_meas(i) = v_store_true_w(sensorNode(i), k+1) + v_noise_std*randn;
    end

    %% ---------------- State prediction (use current fd_cur to build K) --------------
    % Build modal stiffness from current parameter estimate (reset K_dkf each step)
    K_dkf_modal = zeros(modes,modes);
    for z = 1:Nz
        K_dkf_modal = K_dkf_modal + (1 - fd_cur(z)) * Kz_modal_red{z};
    end

    % Continuous-time modal A_cont (correct sign)
    A_cont = [ zeros(modes), eye(modes);
              - (M_modal_red \ K_dkf_modal),  - (M_modal_red \ C_modal_red) ];  % 2m x 2m

    % Continuous-time B_cont maps modal force (F_modal_red) into states
    B_cont = [ zeros(modes); (M_modal_red \ eye(modes)) ];   % 2m x m
    % Discretize (exact) using matrix exponential on augmented matrix
    nA = 2*modes;
    Mbig = [A_cont, B_cont; zeros(modes, nA), zeros(modes)];
    E = expm(Mbig*dt);
    Ad = E(1:nA, 1:nA);
    Bd = E(1:nA, nA+1:end);

    % predict
    x_pred = Ad * x_cur + Bd * F_modal_red;
    P_x_pred = Ad * P_x * Ad' + Q_x;

    %% ---------------- State update (using u_meas and v_meas) -----------------------
    % measurement matrix mapping x=[q; qdot] to [u_meas; v_meas]
    % u = Phi_sen * q, v = Phi_sen * qdot
    Hx = [Phi_sen, zeros(Ns, modes); zeros(Ns, modes), Phi_sen];  % 2Ns x 2m

    Z_uv = [u_meas; v_meas];    % 2Ns x 1
    z_uv_diff = Z_uv - Hx * x_pred;

    Sx = Hx * P_x_pred * Hx' + R_x;
    % regularize Sx to avoid singularity
    Sx = Sx + 1e-12 * eye(size(Sx));
    Kx = (P_x_pred * Hx') / Sx;       % 2m x 2Ns

    x_upd = x_pred + Kx * z_uv_diff;
    P_x = (Istates - Kx * Hx) * P_x_pred;

    x_cur = x_upd;
    x_state(:, k+1) = x_cur;

    %% ---------------- Parameter prediction (random walk) --------------------------
    A_fd = eye(params);
    fd_pred = A_fd * fd_cur;                % no dynamics: random walk
    P_fd_pred = A_fd * P_fd * A_fd' + Q_fd;

    %% ---------------- Parameter update using acceleration measurement -----------
    % Predict acceleration at sensors from state prediction (modal)
    q_pred = x_cur(1:modes);
    qdot_pred = x_cur(modes+1:end);
    qdd_pred = M_modal_red \ (F_modal_red - K_dkf_modal * q_pred - C_modal_red * qdot_pred); % modal accel
    a_pred_sensors = Phi_sen * qdd_pred;   % Ns x 1

    % innovation
    Z_accn = a_meas(:, k+1);
    z_accn_diff = Z_accn - a_pred_sensors;

    % build Jacobian H_fd (Ns x Nz)
    H_fd = zeros(Ns, params);
    for z = 1:params
        % derivative: S * M^{-1} * Kz * u
        % using modal: Phi_sen * (M_modal_red \ (Kz_modal_red{z} * q_pred))
        tmp = M_modal_red \ (Kz_modal_red{z} * q_pred);    % modes x 1
        H_fd(:, z) = Phi_sen * tmp;                        % Ns x 1
    end

    % update parameters using EKF formula
    S_fd = H_fd * P_fd_pred * H_fd' + R_fd;
    S_fd = S_fd + 1e-12 * eye(size(S_fd));  % regularize
    K_fd = (P_fd_pred * H_fd') / S_fd;      % Nz x Ns

    fd_upd = fd_pred + K_fd * z_accn_diff;
    % enforce physical bounds and remove NaNs
    fd_upd = real(fd_upd);                  % remove tiny imaginary parts
    fd_upd(isnan(fd_upd)) = 0;              % safety
    fd_upd = max(0, min(0.99, fd_upd));     % clamp to [0,0.99]

    P_fd = (Iparams - K_fd * H_fd) * P_fd_pred;

    fd_cur = fd_upd;
    fd_est(:, k+1) = fd_cur;

    %% store DKF-predicted scattered fields (for plotting)
    % reconstruct predicted physical displacement at all free DOFs:
    U_modal = Phi_modal * x_cur(1:modes);
    Ufull = zeros(Ndof,1);
    Ufull(freeDOF) = U_modal;
    u_store_dkf_w(:, k+1) = Ufull(w_dofs);

    % predicted accel (reconstruct)
    qdd = qdd_pred;
    Uacc_modal = Phi_modal * qdd;
    Ufull = zeros(Ndof,1);
    Ufull(freeDOF) = Uacc_modal;
    a_store_dkf_w(:, k+1) = Ufull(w_dofs);
end

fprintf('DKF finished.\n');

%% ---------------- PLOTS (essential ones) ----------------
time = tvec;

% 1) True vs estimated damage over time
figure('Name','Damage estimates');
for z = 1:Nz
    subplot(Nz,1,z);
    plot(time, fd(z, :), 'k-', 'LineWidth',1.5); hold on;
    plot(time, fd_est(z, :), 'r--', 'LineWidth',1.5);
    ylabel(sprintf('\\alpha_{%d}', z)); grid on;
    if z==1
        title('True (black) vs Estimated (red) damage factors');
    end
    if z==Nz
        xlabel('Time (s)');
    end
end

% 2) Sensor accelerations: true vs predicted (first sensor example)
figure('Name','Sensor accelerations');
for i = 1:Ns
    subplot(Ns,1,i);
    plot(time, a_store_true_w(sensorNode(i),:), 'k-'); hold on;
    plot(time, a_store_dkf_w(sensorNode(i),:), 'r--');
    ylabel(sprintf('a @ S%d (m/s^2)', i));
    grid on;
    if i==1, title('True (black) vs DKF-predicted (red) sensor accel'); end
    if i==Ns, xlabel('Time (s)'); end
end

% 3) Displacement at a sensor (u) true vs DKF
figure('Name','Sensor displacements');
for i = 1:Ns
    subplot(Ns,1,i);
    plot(time, u_store_true_w(sensorNode(i),:), 'k-'); hold on;
    plot(time, u_store_dkf_w(sensorNode(i),:), 'r--');
    ylabel(sprintf('w @ S%d (m)', i)); grid on;
    if i==1, title('True vs DKF-predicted displacement at sensors'); end
    if i==Ns, xlabel('Time (s)'); end
end

% 4) Modal states (first few)
figure('Name','Modal displacements (first modes)');
for m = 1:min(4,modes)
    subplot(4,1,m);
    plot(time, x_state(m, :));
    ylabel(sprintf('q_%d', m));
    grid on;
    if m==1, title('Estimated modal coordinates'); end
end
xlabel('Time (s)');

% 5) Visual snapshot of deflection at final time (true vs DKF)
figure('Name','Final deflection (true vs dkf)');
% reshape grid for plotting
nxp = nx+1; nyp = ny+1;
X = reshape(coords(:,1), nxp, nyp)';
Y = reshape(coords(:,2), nxp, nyp)';
Wtrue = reshape(u_store_true_w(:, end), nxp, nyp)';
Wdkf  = reshape(u_store_dkf_w(:, end), nxp, nyp)';
subplot(1,2,1);
surf(X, Y, Wtrue); title('True deflection'); shading interp; colorbar; view(45,30);
subplot(1,2,2);
surf(X, Y, Wdkf); title('DKF predicted deflection'); shading interp; colorbar; view(45,30);

fprintf('Plots created.\n');

% Vizualize the real and estimated results ( disp at node center ) 
figure('Name','Displacement plot at sensor');
for i=1:Ns
    subplot(Ns,1,i);
    plot(time,u_store_true_w(sensorNode(i),:),'k-'); hold on;
    plot(time,u_store_dkf_w(sensorNode(i),:),'r--');
    ylabel('Displacment'); grid on;
end  
