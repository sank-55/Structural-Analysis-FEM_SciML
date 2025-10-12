%% Plate_Kirchhoff_FEM_Newmark_KalmanFilter_Enhanced.m
% Compares the KF estimated points on the sensor Measured points then show the defo
% There are some Fault have to look into this thing
% Enhanced version with comprehensive measurement data and improved visualization
% Comparison of FEM Analytical (Black), Noisy Measurements (Red), KF Estimated (Blue)

clear; close all; clc;

%% ---------------- MATERIAL & GEOMETRY ----------------
E = 210e9;         % Young's modulus (Pa)
h = 0.005;         % Plate thickness (m)
rho = 8050;        % Density (kg/m³)
nu = 0.3;          % Poisson's ratio

Lx = 2; Ly = 2;    % Plate dimensions (m)
q0 = -1500;        % Distributed load (N/m²)

% Mesh configuration
nx = 16; ny = 16;  % Balanced mesh for accuracy and speed
dx = Lx/nx; dy = Ly/ny;
Nelem = nx*ny;
Nnode = (nx+1)*(ny+1);
Ndof = 3*Nnode;    % [w, theta_x, theta_y] per node

fprintf('=== PLATE FEM WITH KALMAN FILTER ===\n');
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

% Element connectivity
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

%% ---------------- SHAPE FUNCTIONS & GAUSS INTEGRATION ----------------
gp = [-1/sqrt(3), 1/sqrt(3)];
wgp = [1,1];

syms xi eta real
% Hermite cubic shape functions
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
fprintf('Assembling global matrices...\n');
Kg = sparse(Ndof, Ndof);
Mg = sparse(Ndof, Ndof);
Fg = zeros(Ndof,1);

for e = 1:Nelem
    nodes = c(e,:);
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
    Kg(dofs,dofs) = Kg(dofs,dofs) + Ke;
    Mg(dofs,dofs) = Mg(dofs,dofs) + Me;
    Fg(dofs) = Fg(dofs) + Fe;
end

% Numerical stabilization
Kg = Kg + 1e-10 * speye(size(Kg));

%% ---------------- BOUNDARY CONDITIONS ----------------
tol = 1e-9;
left_nodes = find(abs(coords(:,1)) < tol);
fixedDOF = reshape([3*left_nodes-2; 3*left_nodes-1; 3*left_nodes],1,[]);
freeDOF = setdiff(1:Ndof, fixedDOF);

Kf = Kg(freeDOF, freeDOF);
Mf = Mg(freeDOF, freeDOF);
Ff = Fg(freeDOF);

fprintf('Free DOF: %d, Fixed DOF: %d\n', length(freeDOF), length(fixedDOF));

%% ---------------- MODAL ANALYSIS & RAYLEIGH DAMPING ----------------
nEig = min(10, length(freeDOF)-2);
try
    [PHI_all, W2_all] = eigs(Kf, Mf, nEig, 'sm');
    omega_all = sqrt(abs(diag(W2_all)));
catch
    [PHIs, W2s] = eig(full(Kf), full(Mf));
    omega_all = sqrt(abs(diag(W2s)));
    [omega_all, idxs] = sort(omega_all);
    PHI_all = PHIs(:, idxs(1:min(nEig,length(idxs))));
end

% Rayleigh damping
zeta_target = 0.02;
if length(omega_all) >= 2
    w1 = omega_all(1); w2 = omega_all(2);
    A = [1/(2*w1), w1/2; 1/(2*w2), w2/2];
    bvec = [zeta_target; zeta_target];
    ab = A\bvec;
    alphaR = ab(1); betaR = ab(2);
else
    alphaR = 0.1; betaR = 0.001;
end
Cf = alphaR * Mf + betaR * Kf;
fprintf('Natural frequencies: %.2f, %.2f, %.2f Hz\n', omega_all(1:3)/(2*pi));
fprintf('Rayleigh damping: alpha=%.2e, beta=%.2e\n', alphaR, betaR);

%% ---------------- KALMAN FILTER SETUP WITH EXTENSIVE MEASUREMENTS ----------------
fprintf('\n=== KALMAN FILTER CONFIGURATION ===\n');

% Enhanced measurement parameters
measurement_noise_std = 2e-6;  % Reduced noise for better estimation less is near the true
n_measured_nodes = min(25, floor(Nnode/2));  % Increased measurement points

% Strategic measurement node selection - cover entire plate
all_nodes = 1:Nnode;
free_nodes = setdiff(all_nodes, left_nodes);

% Select measurement nodes in grid pattern for better coverage
x_coords = coords(free_nodes,1);
y_coords = coords(free_nodes,2);

% Create measurement grid
x_meas_grid = linspace(min(x_coords), max(x_coords), 5);
y_meas_grid = linspace(min(y_coords), max(y_coords), 5);
[X_meas, Y_meas] = meshgrid(x_meas_grid, y_meas_grid);
measurement_grid = [X_meas(:), Y_meas(:)];

% Find closest nodes to grid points
measured_node_indices = [];
for i = 1:size(measurement_grid,1)
    distances = sqrt((coords(free_nodes,1)-measurement_grid(i,1)).^2 + ...
                    (coords(free_nodes,2)-measurement_grid(i,2)).^2);
    [~, idx] = min(distances);
    measured_node_indices = [measured_node_indices, free_nodes(idx)];
end
measured_node_indices = unique(measured_node_indices);
measured_node_indices = measured_node_indices(1:min(n_measured_nodes, length(measured_node_indices)));

fprintf('Measurement nodes: %d (%.1f%% coverage)\n', length(measured_node_indices), ...
    100*length(measured_node_indices)/Nnode);

% Measurement matrix H
H = zeros(length(measured_node_indices), length(freeDOF));
valid_measurements = 0;

for i = 1:length(measured_node_indices)
    node_idx = measured_node_indices(i);
    dof_idx = 3*(node_idx-1) + 1;  % w displacement DOF
    col_idx = find(freeDOF == dof_idx, 1);
    if ~isempty(col_idx)
        valid_measurements = valid_measurements + 1;
        H(valid_measurements, col_idx) = 1;
    end
end

H = H(1:valid_measurements, :);
fprintf('Valid measurement DOF: %d\n', valid_measurements);

% Kalman Filter parameters
n_states = length(freeDOF);
Q_kf = 1e-4 * speye(n_states);  % Process noise  MORE IS MORE BIAS TO NOISY DATA AND VICE VERCA
R_kf = (measurement_noise_std^2) * speye(size(H,1));  % Measurement noise

% Initial state and covariance
x_hat = zeros(n_states, 1);
P_kf = 1e-8 * speye(n_states);

% Storage for analysis
kf_errors = [];
kf_gain_history = [];
innovation_history = [];

%% ---------------- NEWMARK TIME INTEGRATION ----------------
beta_nm = 1/4; gamma_nm = 1/2;
dt = 0.001;      % Time step
Ttotal = 1.0;    % Total simulation time
tvec = 0:dt:Ttotal;
nSteps = length(tvec)-1;

Kred = Kf; Mred = Mf; Cred = Cf;
F_static = Ff;

% Ramped loading
tramp = 0.1;
load_factor = @(tt) min(1, tt/tramp);

% Newmark constants
a0 = 1/(beta_nm*dt^2);
a1 = gamma_nm/(beta_nm*dt);
a2 = 1/(beta_nm*dt);
a3 = 1/(2*beta_nm)-1;
a4 = gamma_nm/beta_nm - 1;
a5 = dt*(gamma_nm/(2*beta_nm) - 1);

Keff = Kred + a0*Mred + a1*Cred;
Keff = Keff + 1e-12 * speye(size(Keff));  % Regularization

% LU decomposition for efficient solving
[L, U, P, Q] = lu(Keff);

% Initial conditions
nFree = size(Kred,1);
u = zeros(nFree,1);
v = zeros(nFree,1);
a = Mred \ ( (F_static * load_factor(0)) - Cred*v - Kred*u );

% Initialize Kalman Filter with first measurement
initial_measurement = H * u;
initial_noise = measurement_noise_std * randn(size(initial_measurement));
x_hat = H' * ((H * H' + 1e-14 * eye(size(H,1))) \ (initial_measurement + initial_noise));

% Storage arrays
Ufree_hist = zeros(nFree, nSteps+1);
w_time = zeros(Nnode, nSteps+1);
measurements_hist = zeros(size(H,1), nSteps+1);
kf_w_time = zeros(Nnode, nSteps+1);

% Store initial conditions
Ufree_hist(:,1) = u;
Ufull = zeros(Ndof,1);
Ufull(freeDOF) = u;
w_time(:,1) = Ufull(1:3:end);
measurements_hist(:,1) = initial_measurement + initial_noise;

kf_full = zeros(Ndof,1);
kf_full(freeDOF) = x_hat;
kf_w_time(:,1) = kf_full(1:3:end);

fprintf('\n=== TIME INTEGRATION STARTED ===\n');
fprintf('Time steps: %d, dt: %.4f s, Total time: %.2f s\n', nSteps, dt, Ttotal);

tic;
for step = 1:nSteps
    tn1 = tvec(step+1);
    
    %% TRUE FEM SOLUTION
    Fext_n1 = F_static * load_factor(tn1);
    Meff_rhs = Fext_n1 + Mred*(a0*u + a2*v + a3*a) + Cred*(a1*u + a4*v + a5*a);
    
    u_new = Q * (U \ (L \ (P * Meff_rhs)));
    a_new = a0*(u_new - u) - a2*v - a3*a;
    v_new = v + dt*((1-gamma_nm)*a + gamma_nm*a_new);
    
    % Store true solution
    Ufree_hist(:, step+1) = u_new;
    Ufull = zeros(Ndof,1);
    Ufull(freeDOF) = u_new;
    w_time(:, step+1) = Ufull(1:3:end);
    
    %% GENERATE NOISY MEASUREMENTS
    true_measurement = H * u_new;
    noisy_measurement = true_measurement + measurement_noise_std * randn(size(true_measurement));
    measurements_hist(:, step+1) = noisy_measurement;
    
    %% KALMAN FILTER
    % Prediction
    A_kf = speye(n_states);
    x_pred = A_kf * x_hat;
    P_pred = A_kf * P_kf * A_kf' + Q_kf;
    P_pred = (P_pred + P_pred')/2 + 1e-14 * speye(size(P_pred));
    
    % Update
    S = H * P_pred * H' + R_kf;
    S = (S + S')/2 + 1e-14 * speye(size(S));
    
    try
        K_gain = P_pred * H' / S;
    catch
        K_gain = P_pred * H' * pinv(S);
    end
    
    innovation = noisy_measurement - H * x_pred;
    x_hat = x_pred + K_gain * innovation;
    
    % Joseph form covariance update
    I = speye(n_states);
    P_kf = (I - K_gain * H) * P_pred * (I - K_gain * H)' + K_gain * R_kf * K_gain';
    P_kf = (P_kf + P_kf')/2;
    
    % Store KF results
    kf_full = zeros(Ndof,1);
    kf_full(freeDOF) = x_hat;
    kf_w_time(:, step+1) = kf_full(1:3:end);
    
    % Store diagnostics
    kf_error = norm(u_new - x_hat) / (norm(u_new) + 1e-12);
    kf_errors = [kf_errors, kf_error];
    kf_gain_history = [kf_gain_history, mean(diag(K_gain * K_gain'))];
    innovation_history = [innovation_history, norm(innovation)];
    
    % Update for next step
    u = u_new; v = v_new; a = a_new;
    
    if mod(step, max(1,round(nSteps/20)))==0
        fprintf('  Progress: %3d%%, t=%.3f s, KF Error: %.2e\n', ...
            round(100*step/nSteps), tn1, kf_error);
    end
end
computation_time = toc;
fprintf('Time integration completed in %.2f seconds\n', computation_time);

%% ---------------- COMPREHENSIVE RESULTS ANALYSIS ----------------
fprintf('\n=== RESULTS ANALYSIS ===\n');

% Key node indices
[~, centerIdx] = min((coords(:,1)-Lx/2).^2 + (coords(:,2)-Ly/2).^2);
[~, quarterIdx] = min((coords(:,1)-Lx/4).^2 + (coords(:,2)-Ly/4).^2);
[~, threeQuarterIdx] = min((coords(:,1)-3*Lx/4).^2 + (coords(:,2)-3*Ly/4).^2);

key_nodes = [centerIdx, quarterIdx, threeQuarterIdx];
node_names = {'Center', 'Quarter', 'Three-Quarter'};

% Extract data for key nodes
node_data = struct();
for i = 1:length(key_nodes)
    node_idx = key_nodes(i);
    node_data(i).name = node_names{i};
    node_data(i).fem = w_time(node_idx, :);
    node_data(i).kf = kf_w_time(node_idx, :);
    
    % Find if this node is measured
    is_measured = ismember(node_idx, measured_node_indices);
    if is_measured
        [~, meas_idx] = ismember(node_idx, measured_node_indices);
        node_data(i).measured = measurements_hist(meas_idx, :);
    else
        node_data(i).measured = nan(size(node_data(i).fem));
    end
    
    % Calculate errors
    node_data(i).kf_error = abs(node_data(i).fem - node_data(i).kf);
    if is_measured
        node_data(i).meas_error = abs(node_data(i).fem - node_data(i).measured);
    end
end

% Calculate RMSE values
rmse_results = [];
for i = 1:length(key_nodes)
    valid_kf = ~isnan(node_data(i).kf) & ~isinf(node_data(i).kf);
    rmse_kf = sqrt(mean(node_data(i).kf_error(valid_kf).^2));
    
    if ~all(isnan(node_data(i).measured))
        valid_meas = ~isnan(node_data(i).measured) & ~isinf(node_data(i).measured);
        rmse_meas = sqrt(mean(node_data(i).meas_error(valid_meas).^2));
    else
        rmse_meas = nan;
    end
    
    rmse_results = [rmse_results; rmse_meas, rmse_kf];
    fprintf('%s Node - RMSE Measurement: %.3e m, RMSE KF: %.3e m\n', ...
        node_data(i).name, rmse_meas, rmse_kf);
end

%% ---------------- ENHANCED VISUALIZATION ----------------
fprintf('\nGenerating comprehensive plots...\n');

% FIGURE 1: MULTI-NODE COMPARISON
figure('Position', [50, 100, 1400, 900]);

for i = 1:length(key_nodes)
    subplot(2, 2, i);
    
    % Main comparison plot
    plot(tvec, node_data(i).fem, 'k-', 'LineWidth', 3, 'DisplayName', 'FEM Analytical');
    hold on;
    
    if ~all(isnan(node_data(i).measured))
        plot(tvec, node_data(i).measured, 'ro', 'MarkerSize', 4, ...
            'MarkerFaceColor', 'r', 'DisplayName', 'Noisy Measurements');
    end
    
    plot(tvec, node_data(i).kf, 'b-', 'LineWidth', 2, 'DisplayName', 'KF Estimated');
    
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Displacement w (m)', 'FontSize', 11);
    title(sprintf('%s Node: FEM vs Measurements vs KF', node_data(i).name), 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 10);
end

% FIGURE 1d: Error comparison
subplot(2, 2, 4);
error_data = [];
labels = {};
for i = 1:length(key_nodes)
    if ~all(isnan(node_data(i).measured))
        error_data = [error_data, rmse_results(i,1)];
        labels{end+1} = [node_data(i).name ' Meas'];
    end
    error_data = [error_data, rmse_results(i,2)];
    labels{end+1} = [node_data(i).name ' KF'];
end

bar(error_data, 'FaceColor', [0.6 0.6 0.6]);
set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 45);
ylabel('RMS Error (m)');
title('RMS Error Comparison Across Key Nodes');
grid on;

sgtitle('Comprehensive Node-wise Comparison: FEM (Black) vs Measurements (Red) vs KF (Blue)', ...
    'FontSize', 14, 'FontWeight', 'bold');

% FIGURE 2: KALMAN FILTER PERFORMANCE METRICS
figure('Position', [100, 50, 1200, 800]);

subplot(2,3,1);
plot(tvec(2:end), kf_errors, 'k-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Normalized State Error');
title('Kalman Filter Error Convergence');
grid on;

subplot(2,3,2);
plot(tvec(2:end), kf_gain_history, 'm-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Average Kalman Gain');
title('Kalman Gain Evolution');
grid on;

subplot(2,3,3);
plot(tvec(2:end), innovation_history, 'g-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Innovation Norm');
title('Innovation Sequence');
grid on;

subplot(2,3,4);
% Spatial error at final time
final_error = abs(w_time(:,end) - kf_w_time(:,end));
scatter3(coords(:,1), coords(:,2), final_error, 40, final_error, 'filled');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Error (m)');
title('Spatial Error Distribution (Final Time)');
colorbar; grid on;

subplot(2,3,5);
% Measurement locations
plot(coords(:,1), coords(:,2), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'k');
hold on;
plot(coords(measured_node_indices,1), coords(measured_node_indices,2), 'ro', ...
     'MarkerSize', 6, 'MarkerFaceColor', 'r');
for i = 1:length(key_nodes)
    plot(coords(key_nodes(i),1), coords(key_nodes(i),2), 'bs', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'b');
end
xlabel('X (m)'); ylabel('Y (m)');
title('Sensor Network (Red) and Key Nodes (Blue)');
legend('All Nodes', 'Measurement Nodes', 'Key Analysis Nodes', 'Location', 'best');
axis equal; grid on;

subplot(2,3,6);
% Performance improvement
improvement = (rmse_results(:,1) - rmse_results(:,2)) ./ rmse_results(:,1) * 100;
bar(improvement, 'FaceColor', [0.2 0.6 0.8]);
set(gca, 'XTickLabel', node_names);
ylabel('Error Reduction (%)');
title('KF Performance Improvement Over Measurements');
grid on;

sgtitle('Kalman Filter Performance Analysis', 'FontSize', 14, 'FontWeight', 'bold');

% FIGURE 3: ANIMATION FRAME COMPARISON
figure('Position', [150, 80, 1000, 800]);

% Final time step comparison
final_fem = w_time(:,end);
final_kf = kf_w_time(:,end);

% Create interpolated surfaces
xq = linspace(0, Lx, 50);
yq = linspace(0, Ly, 50);
[Xq, Yq] = meshgrid(xq, yq);

F_fem = scatteredInterpolant(coords(:,1), coords(:,2), final_fem, 'natural', 'none');
F_kf = scatteredInterpolant(coords(:,1), coords(:,2), final_kf, 'natural', 'none');

Z_fem = F_fem(Xq, Yq);
Z_kf = F_kf(Xq, Yq);

subplot(2,2,1);
surf(Xq, Yq, Z_fem, 'EdgeColor', 'none');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('w (m)');
title('FEM Analytical Solution (Final Time)');
colorbar; view(45,30); axis equal;

subplot(2,2,2);
surf(Xq, Yq, Z_kf, 'EdgeColor', 'none');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('w (m)');
title('KF Estimated Solution (Final Time)');
colorbar; view(45,30); axis equal;

subplot(2,2,3);
error_surface = abs(Z_fem - Z_kf);
surf(Xq, Yq, error_surface, 'EdgeColor', 'none');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Error (m)');
title('Estimation Error (Final Time)');
colorbar; view(45,30); axis equal;

subplot(2,2,4);
% Measurement points on error surface
scatter3(coords(measured_node_indices,1), coords(measured_node_indices,2), ...
         final_error(measured_node_indices), 50, 'r', 'filled');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Error (m)');
title('Measurement Points Error Distribution');
grid on; view(45,30);

sgtitle('Final Time Step Spatial Comparison', 'FontSize', 14, 'FontWeight', 'bold');

%% ---------------- PERFORMANCE SUMMARY ----------------
fprintf('\n=== KALMAN FILTER PERFORMANCE SUMMARY ===\n');
fprintf('Configuration:\n');
fprintf('  Measurement nodes: %d/%d (%.1f%% coverage)\n', ...
    length(measured_node_indices), Nnode, 100*length(measured_node_indices)/Nnode);
fprintf('  Measurement noise: %.2e m\n', measurement_noise_std);
fprintf('  Computation time: %.2f seconds\n', computation_time);

fprintf('\nConvergence Analysis:\n');
if ~isempty(kf_errors)
    initial_error = kf_errors(1);
    final_error = kf_errors(end);
    convergence_ratio = final_error / initial_error;
    
    fprintf('  Initial KF error: %.3e\n', initial_error);
    fprintf('  Final KF error: %.3e\n', final_error);
    fprintf('  Convergence ratio: %.3f\n', convergence_ratio);
    
    if convergence_ratio < 0.01
        fprintf(' Excellent convergence (%.1f%% of initial error)\n', convergence_ratio*100);
    elseif convergence_ratio < 0.1
        fprintf(' Very good convergence (%.1f%% of initial error)\n', convergence_ratio*100);
    elseif convergence_ratio < 0.3
        fprintf(' Good convergence (%.1f%% of initial error)\n', convergence_ratio*100);
    else
        fprintf(' Moderate convergence (%.1f%% of initial error)\n', convergence_ratio*100);
    end
end

fprintf('\nOverall Performance:\n');
avg_rmse_kf = mean(rmse_results(:,2));
fprintf('  Average KF RMSE across key nodes: %.3e m\n', avg_rmse_kf);
fprintf('  Maximum KF error reduction: %.1f%%\n', max(improvement));

fprintf('\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n');

% Save key results
results = struct();
results.tvec = tvec;
results.coords = coords;
results.w_fem = w_time;
results.w_kf = kf_w_time;
results.measurements = measurements_hist;
results.measured_nodes = measured_node_indices;
results.kf_errors = kf_errors;
results.rmse_results = rmse_results;

fprintf('Results structure saved for further analysis.\n');
