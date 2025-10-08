%% plate_cantilever_kf_smooth.m
% Minimal cantilever plate + smooth forcing + Kalman filter
clear; close all; clc; rng(0);

%% ---------------- User / mesh parameters ----------------
E = 210e9; nu = 0.3; h = 0.01; rho = 7800;
Lx = 1.0; Ly = 0.5;
nx = 6; ny = 3;                 % keep small for demo
dx = Lx/nx; dy = Ly/ny;
Nnode = (nx+1)*(ny+1);
Ndof = 3 * Nnode;

% element scale (simple demo, not full Kirchhoff matrices)
k_elem_scale = 5e7;    % larger -> stiffer -> smoother dynamics (reduce overshoot)
m_elem_scale = rho*h*dx*dy/4;

%% ---------------- Mesh & connectivity ----------------
coords = zeros(Nnode,2); node=1;
for j=0:ny, for i=0:nx, coords(node,:)=[i*dx,j*dy]; node=node+1; end, end
c = zeros(nx*ny,4); e=0;
for j=1:ny
  for i=1:nx
    n1=(j-1)*(nx+1)+i; n2=n1+1; n3=n2+(nx+1); n4=n1+(nx+1);
    e=e+1; c(e,:)=[n1 n2 n3 n4];
  end
end
Nelem = size(c,1);

%% ---------------- Assemble simple Kg, Mg ----------------
Kg = zeros(Ndof); Mg = zeros(Ndof);
for e = 1:Nelem
    nodes = c(e,:);
    dofs_e = reshape([3*nodes-2;3*nodes-1;3*nodes],1,[]);
    Ke = k_elem_scale * eye(12);
    Me = m_elem_scale * eye(12);
    Kg(dofs_e,dofs_e) = Kg(dofs_e,dofs_e) + Ke;
    Mg(dofs_e,dofs_e) = Mg(dofs_e,dofs_e) + Me;
end

%% ---------------- Cantilever BC at x=0 ----------------
tol = 1e-12;
fixed_nodes = find(abs(coords(:,1)) < tol);
fixedDOF = reshape([3*fixed_nodes-2;3*fixed_nodes-1;3*fixed_nodes],1,[]);
freeDOF = setdiff(1:Ndof, fixedDOF);
Kf = Kg(freeDOF, freeDOF); Mf = Mg(freeDOF, freeDOF); nFree = length(freeDOF);

%% ---------------- Forcing on right edge ----------------
right_nodes = find(abs(coords(:,1)-Lx) < tol);
F_full = zeros(Ndof,1);
F_full(3*right_nodes - 2) = 1.0;   % spatial shape
Ff_static = F_full(freeDOF);

%% ---------------- Damping (increase slightly for smoothness) ----------------
alphaR = 5e-4;    % increased
betaR  = 5e-6;
Cfull = alphaR * Mg + betaR * Kg;
Cf = Cfull(freeDOF, freeDOF);

%% ---------------- Time integration (Newmark) ----------------
dt = 0.002; Ttotal = 0.6; t = 0:dt:Ttotal; nSteps = length(t);
betaN = 1/4; gammaN = 1/2;
U_true = zeros(nFree, nSteps); V_true = zeros(nFree, nSteps); A_true = zeros(nFree, nSteps);
U_true(:,1)=0; V_true(:,1)=0;
A_true(:,1) = Mf \ (Ff_static*0 - Cf*V_true(:,1) - Kf*U_true(:,1));
Keff = Mf + betaN*dt^2*Kf + gammaN*dt*Cf;

% process noise & measurement noise reduced for smooth plots
proc_acc_std = 1e-5;    % smaller
sigma_meas = 2e-5;      % smaller

% Smooth forcing: sinusoid * exponential ramp (tau small)
f_freq = 5.0; amp = 2.0; tau = 0.03;
force_time = @(tt) amp * (1 - exp(-tt/tau)) .* sin(2*pi*f_freq*tt);  % smooth ramped sine

for k = 1:nSteps-1
    tk = t(k);
    Fext = Ff_static * force_time(tk);
    w_proc = proc_acc_std * randn(nFree,1);
    Fext_noisy = Fext + Mf * w_proc;
    RHS = Fext_noisy + Mf*( (1/dt^2)*U_true(:,k) + (1/dt)*V_true(:,k) + (1/2 - betaN)*A_true(:,k) ) ...
                  + Cf*( (gammaN/(betaN*dt))*U_true(:,k) + ((gammaN/betaN)-1)*V_true(:,k) + dt*((gammaN/(2*betaN))-1)*A_true(:,k) );
    Unew = Keff \ RHS;
    Anew = (1/(betaN*dt^2))*(Unew - U_true(:,k) - dt*V_true(:,k)) - (1/(betaN*dt))*V_true(:,k) - ((1/(2*betaN))-1)*A_true(:,k);
    Vnew = V_true(:,k) + dt*((1-gammaN)*A_true(:,k) + gammaN*Anew);
    U_true(:,k+1)=Unew; V_true(:,k+1)=Vnew; A_true(:,k+1)=Anew;
end

%% ---------------- Sensor: mid node on right edge ----------------
right_coords = coords(right_nodes,:);
[~, idx_mid] = min(sum((right_coords - repmat([Lx, Ly/2], size(right_coords,1),1)).^2,2));
sensor_node = right_nodes(idx_mid);
sensor_full_dof = 3*sensor_node - 2;
sensor_free_index = find(freeDOF == sensor_full_dof, 1);
if isempty(sensor_free_index), error('Sensor on fixed DOF'); end

% measurements (smooth + small noise)
y_meas = U_true(sensor_free_index, :) + sigma_meas * randn(1,nSteps);

%% ---------------- State-space (discrete Euler) for KF ----------------
I_n = speye(nFree); Minv = Mf \ I_n;
Ac = [sparse(nFree,nFree), I_n; -Minv*Kf, -Minv*Cf];
Bc = [sparse(nFree,nFree); Minv];
nState = 2*nFree;
Fd = eye(nState) + full(Ac)*dt;
Gd = full(Bc)*dt;

Qc = (proc_acc_std^2) * eye(nFree);
Qd = Gd * Qc * Gd';

H = zeros(1, nState); H(1, sensor_free_index) = 1;   % measure that displacement DOF
Rk = sigma_meas^2;

%% ---------------- Kalman Filter ----------------
Xkf = zeros(nState, nSteps);
Pk = repmat(1e-6 * eye(nState), [1,1,nSteps]);
Xkf(:,1) = 0; Pk(:,:,1) = 1e-4 * eye(nState);

for k = 1:nSteps-1
    fk = Ff_static * force_time(t(k));          % nFree x 1
    Xbar = Fd * Xkf(:,k) + Gd * fk;
    Pbar = Fd * Pk(:,:,k) * Fd' + Qd;
    yk = y_meas(k+1);
    S = H * Pbar * H' + Rk;
    Kk = (Pbar * H') / S;
    innov = yk - H * Xbar;
    Xpost = Xbar + Kk * innov;
    Ppost = Pbar - Kk * H * Pbar;
    Xkf(:,k+1) = Xpost; Pk(:,:,k+1) = Ppost;
end

%% ---------------- Final single plot (smooth true, measurements, KF) ----------------
true_ts = U_true(sensor_free_index, :);
kf_ts   = Xkf(sensor_free_index, :);
meas_ts = y_meas;

figure('Color','w','Position',[200 200 800 350]);
plot(t, true_ts, 'k-', 'LineWidth', 1.6); hold on;
plot(t, meas_ts, 'r.', 'MarkerSize', 8);
plot(t, kf_ts, 'b--', 'LineWidth', 1.4);
legend('True','Measurement','KF estimate','Location','best');
xlabel('Time (s)'); ylabel('Transverse deflection at sensor (m)');
title('Sensor: True (smooth) vs Measurement vs KF estimate');
grid on; box on;
