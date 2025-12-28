%% DUAL EXTENDED KALMAN FILTER FOR 3-DOF SYSTEM WITH MULTIPLE DAMAGES
% Simultaneous tracking of progressive damage in all three springs

%% ---------------------------------
%     This snippet is able to predict the states and damage parameters very nicely one of the greatest dkf 
%% -------------------------------------------------------

clear; close all; clc;

%% ========================================================================
% SYSTEM DEFINITION - 3-DOF MASS-SPRING-DAMPER
% ========================================================================
fprintf('=== 3-DOF Mass-Spring-Damper System with Multiple Damages ===\n');

% True system parameters
m1 = 10; m2 = 8; m3 = 6;       % Masses (kg)
k1 = 1000; k2 = 800; k3 = 600; % Initial spring stiffness (N/m)
c1 = 8; c2 = 6; c3 = 4;        % Damping coefficients (Ns/m)

% Mass matrix
M = diag([m1, m2, m3]);

% Initial stiffness matrix
K_initial = [k1 + k2, -k2, 0;
             -k2, k2 + k3, -k3;
             0, -k3, k3];

% Initial damping matrix (proportional damping)
C_initial = [c1 + c2, -c2, 0;
             -c2, c2 + c3, -c3;
             0, -c3, c3];

fprintf('Initial parameters:\n');
fprintf('k1 = %.0f N/m, k2 = %.0f N/m, k3 = %.0f N/m\n', k1, k2, k3);
fprintf('c1 = %.1f Ns/m, c2 = %.1f Ns/m, c3 = %.1f Ns/m\n', c1, c2, c3);

%% ========================================================================
% TRUE DAMAGE SCENARIO - ALL SPRINGS DAMAGED PROGRESSIVELY
% ========================================================================
% Time parameters
dt = 0.001;          % Smaller time step for better tracking
T = 15;              % Total simulation time (s)
N = round(T/dt);     % Number of time steps
time = linspace(0, T, N+1);

% Define progressive damage for each spring
% Damage factor: 0 = no damage, 1 = complete failure

% Spring 1: Early damage, quick progression
damage_k1 = zeros(1, N+1);
for i = 1:N+1
    if time(i) <= 3
        damage_k1(i) = 0.4 * (time(i)/3);  % 40% damage by 3s
    elseif time(i) <= 7
        damage_k1(i) = 0.4 + 0.2 * ((time(i)-3)/4); % Additional 20% by 7s
    else
        damage_k1(i) = 0.6;  % Total 60% damage
    end
end

% Spring 2: Later damage, slower progression
damage_k2 = zeros(1, N+1);
for i = 1:N+1
    if time(i) <= 5
        damage_k2(i) = 0;  % No damage initially
    elseif time(i) <= 10
        damage_k2(i) = 0.5 * ((time(i)-5)/5);  % 50% damage by 10s
    else
        damage_k2(i) = 0.5;  % Total 50% damage
    end
end

% Spring 3: Gradual continuous damage
damage_k3 = zeros(1, N+1);
for i = 1:N+1
    if time(i) <= 8
        damage_k3(i) = 0.3 * (time(i)/8);  % 30% damage by 8s
    elseif time(i) <= 12
        damage_k3(i) = 0.3 + 0.1 * ((time(i)-8)/4); % Additional 10% by 12s
    else
        damage_k3(i) = 0.4;  % Total 40% damage
    end
end

% True stiffness values over time
k1_true = k1 * (1 - damage_k1);
k2_true = k2 * (1 - damage_k2);
k3_true = k3 * (1 - damage_k3);

% Damping remains constant (or could also be damaged)
c1_true = c1 * ones(1, N+1);
c2_true = c2 * ones(1, N+1);
c3_true = c3 * ones(1, N+1);

fprintf('\nDamage scenario:\n');
fprintf('Spring 1: Progressive damage up to %.0f%%\n', max(damage_k1)*100);
fprintf('Spring 2: Progressive damage up to %.0f%%\n', max(damage_k2)*100);
fprintf('Spring 3: Progressive damage up to %.0f%%\n', max(damage_k3)*100);

%% ========================================================================
% SIMULATE TRUE SYSTEM RESPONSE WITH PROGRESSIVE DAMAGE
% ========================================================================
fprintf('\nSimulating true system response with progressive damage...\n');

% State vector: [x1; v1; x2; v2; x3; v3]
state_true = zeros(6, N+1);

% External force (colored noise - more realistic)
rng(42);  % For reproducibility
force_freq = 2;  % Hz
F_ext = zeros(3, N+1);
for i = 1:N+1
    % Combine multiple frequency components
    F_ext(:,i) = 50 * (sin(2*pi*force_freq*time(i)) + ...
                       0.3*sin(2*pi*2*force_freq*time(i)) + ...
                       0.1*randn(3,1));
end

% Time integration using Newmark-beta (more stable for varying stiffness)
beta_nm = 0.25;
gamma_nm = 0.5;

% Initialize state
state_true(:,1) = [0.01; 0; 0.005; 0; 0.002; 0];  % Small initial displacements

for i = 1:N
    % Current stiffness values
    k1_curr = k1_true(i);
    k2_curr = k2_true(i);
    k3_curr = k3_true(i);
    
    % Current stiffness matrix
    K_curr = [k1_curr + k2_curr, -k2_curr, 0;
              -k2_curr, k2_curr + k3_curr, -k3_curr;
              0, -k3_curr, k3_curr];
    
    % Current state
    x_curr = state_true(1:2:6, i);
    v_curr = state_true(2:2:6, i);
    
    % Compute acceleration
    a_curr = M \ (F_ext(:,i) - C_initial*v_curr - K_curr*x_curr);
    
    % Predict displacement and velocity
    x_pred = x_curr + dt*v_curr + (0.5-beta_nm)*dt^2*a_curr;
    v_pred = v_curr + (1-gamma_nm)*dt*a_curr;
    
    % Next stiffness values
    k1_next = k1_true(i+1);
    k2_next = k2_true(i+1);
    k3_next = k3_true(i+1);
    
    K_next = [k1_next + k2_next, -k2_next, 0;
              -k2_next, k2_next + k3_next, -k3_next;
              0, -k3_next, k3_next];
    
    % Solve for next acceleration
    a_next = M \ (F_ext(:,i+1) - C_initial*v_pred - K_next*x_pred);
    
    % Update velocity and displacement
    v_next = v_pred + gamma_nm*dt*a_next;
    x_next = x_pred + beta_nm*dt^2*a_next;
    
    % Store
    state_true(:, i+1) = [x_next(1); v_next(1); x_next(2); v_next(2); x_next(3); v_next(3)];
end

%% ========================================================================
% GENERATE NOISY MEASUREMENTS
% ========================================================================
fprintf('Generating noisy measurements...\n');

% Measurement noise levels
SNR_position = 0.0002;   % 2% noise for positions
SNR_velocity = 0.0003;   % 3% noise for velocities
SNR_acceleration = 0.0004; % 4% noise for accelerations

% Extract true states
x1_true = state_true(1,:); v1_true = state_true(2,:);
x2_true = state_true(3,:); v2_true = state_true(4,:);
x3_true = state_true(5,:); v3_true = state_true(6,:);

% Calculate true accelerations
a1_true = zeros(1, N+1);
a2_true = zeros(1, N+1);
a3_true = zeros(1, N+1);
for i = 1:N+1
    k1_curr = k1_true(i);
    k2_curr = k2_true(i);
    k3_curr = k3_true(i);
    K_curr = [k1_curr + k2_curr, -k2_curr, 0;
              -k2_curr, k2_curr + k3_curr, -k3_curr;
              0, -k3_curr, k3_curr];
    
    x_curr = [x1_true(i); x2_true(i); x3_true(i)];
    v_curr = [v1_true(i); v2_true(i); v3_true(i)];
    a_curr = M \ (F_ext(:,i) - C_initial*v_curr - K_curr*x_curr);
    
    a1_true(i) = a_curr(1);
    a2_true(i) = a_curr(2);
    a3_true(i) = a_curr(3);
end

% Add noise to measurements
noise_std_x1 = sqrt(mean(x1_true.^2)) * SNR_position;
noise_std_x2 = sqrt(mean(x2_true.^2)) * SNR_position;
noise_std_x3 = sqrt(mean(x3_true.^2)) * SNR_position;

noise_std_v1 = sqrt(mean(v1_true.^2)) * SNR_velocity;
noise_std_v2 = sqrt(mean(v2_true.^2)) * SNR_velocity;
noise_std_v3 = sqrt(mean(v3_true.^2)) * SNR_velocity;

noise_std_a1 = sqrt(mean(a1_true.^2)) * SNR_acceleration;
noise_std_a2 = sqrt(mean(a2_true.^2)) * SNR_acceleration;
noise_std_a3 = sqrt(mean(a3_true.^2)) * SNR_acceleration;

% Noisy measurements
x1_meas = x1_true + noise_std_x1 * randn(size(x1_true));
x2_meas = x2_true + noise_std_x2 * randn(size(x2_true));
x3_meas = x3_true + noise_std_x3 * randn(size(x3_true));

v1_meas = v1_true + noise_std_v1 * randn(size(v1_true));
v2_meas = v2_true + noise_std_v2 * randn(size(v2_true));
v3_meas = v3_true + noise_std_v3 * randn(size(v3_true));

a1_meas = a1_true + noise_std_a1 * randn(size(a1_true));
a2_meas = a2_true + noise_std_a2 * randn(size(a2_true));
a3_meas = a3_true + noise_std_a3 * randn(size(a3_true));

fprintf('Measurement noise added:\n');
fprintf('  Positions: %.1f%%\n', SNR_position*100);
fprintf('  Velocities: %.1f%%\n', SNR_velocity*100);
fprintf('  Accelerations: %.1f%%\n', SNR_acceleration*100);

%% ========================================================================
% DUAL EXTENDED KALMAN FILTER - MULTIPLE DAMAGE TRACKING
% ========================================================================
fprintf('\n=== DUAL EXTENDED KALMAN FILTER - MULTIPLE DAMAGE TRACKING ===\n');

% Parameters to estimate: [k1, k2, k3, c1, c2, c3]
n_params = 6;
n_states = 6;  % [x1, v1, x2, v2, x3, v3]

% Initial parameter estimates (slightly different from true)
params_est = zeros(n_params, N+1);
params_est(:,1) = [k1*1.1; k2*0.9; k3*1.05; c1*0.8; c2*1.2; c3*0.9];  % 10% error initially

% State estimates
state_est = zeros(n_states, N+1);
state_est(:,1) = [x1_meas(1); v1_meas(1); x2_meas(1); v2_meas(1); x3_meas(1); v3_meas(1)];

% Covariance matrices
P_state = eye(n_states) * 1e-3;     % State covariance
P_param = eye(n_params) * 1e-2;     % Parameter covariance

% Process noise covariance (tuned for better tracking)
Q_state = diag([1e-8, 1e-6, 1e-8, 1e-6, 1e-5, 1e-6]);
Q_param = diag([1e-4, 1e-4, 1e-4, 1e-6, 1e-6, 1e-6]);  % More process noise for stiffness

% Measurement noise covariance
R_state = diag([noise_std_x1^2, noise_std_v1^2, noise_std_x2^2, noise_std_v2^2, ...
                noise_std_x3^2, noise_std_v3^2]);

% Measurement matrix (we measure all states directly)
H_state = eye(n_states);

% Store stiffness estimates separately for plotting
k1_est = zeros(1, N+1); k1_est(1) = params_est(1,1);
k2_est = zeros(1, N+1); k2_est(1) = params_est(2,1);
k3_est = zeros(1, N+1); k3_est(1) = params_est(3,1);

% Progress indicator
progress_step = round(N/20);

% Forgetting factor for parameter tracking (helps track time-varying parameters)
lambda = 0.999;  % Forgetting factor (close to 1)

fprintf('Running DEKF for multiple damage tracking...\n');
for k = 2:N+1
    if mod(k, progress_step) == 0
        fprintf('  Progress: %d%%\n', round(100*k/(N+1)));
    end
    
    % ====================================================================
    % STATE PREDICTION STEP
    % ====================================================================
    % Current parameter estimates
    k1_est_curr = params_est(1, k-1);
    k2_est_curr = params_est(2, k-1);
    k3_est_curr = params_est(3, k-1);
    c1_est = params_est(4, k-1);
    c2_est = params_est(5, k-1);
    c3_est = params_est(6, k-1);
    
    % Build matrices with current parameter estimates
    K_est = [k1_est_curr + k2_est_curr, -k2_est_curr, 0;
             -k2_est_curr, k2_est_curr + k3_est_curr, -k3_est_curr;
             0, -k3_est_curr, k3_est_curr];
    
    C_est = [c1_est + c2_est, -c2_est, 0;
             -c2_est, c2_est + c3_est, -c3_est;
             0, -c3_est, c3_est];
    
    % Continuous state matrix
    A_cont = [zeros(3), eye(3);
              -M\K_est, -M\C_est];
    
    % Input matrix
    B_cont = [zeros(3); inv(M)];
    
    % Discretize (matrix exponential for linear system)
    A_disc = expm(A_cont * dt);
    B_disc = (A_cont\(A_disc - eye(6))) * B_cont;  % More accurate discretization
    
    % State prediction
    state_pred = A_disc * state_est(:, k-1) + B_disc * F_ext(:, k-1);
    
    % State covariance prediction with forgetting factor
    P_state_pred = A_disc * P_state * A_disc' / lambda + Q_state;
    
    % ====================================================================
    % STATE UPDATE STEP
    % ====================================================================
    % Measurement vector (positions and velocities)
    Z_state = [x1_meas(k); v1_meas(k); x2_meas(k); v2_meas(k); x3_meas(k); v3_meas(k)];
    
    % Innovation
    innovation_state = Z_state - H_state * state_pred;
    
    % Innovation covariance
    S_state = H_state * P_state_pred * H_state' + R_state;
    
    % Kalman gain for states
    K_state = P_state_pred * H_state' / S_state;
    
    % State update
    state_est(:, k) = state_pred + K_state * innovation_state;
    
    % State covariance update
    P_state = (eye(n_states) - K_state * H_state) * P_state_pred;
    
    % ====================================================================
    % PARAMETER PREDICTION STEP
    % ====================================================================
    % Parameters evolve as random walk with forgetting factor
    params_pred = params_est(:, k-1);
    P_param_pred = P_param / lambda + Q_param;
    
    % ====================================================================
    % PARAMETER UPDATE STEP - USING ACCELERATION MEASUREMENTS
    % ====================================================================
    % Use acceleration measurements for parameter updates (more sensitive to stiffness)
    Z_accel = [a1_meas(k); a2_meas(k); a3_meas(k)];
    
    % Current state estimates for sensitivity calculation
    x1_curr = state_est(1, k);
    x2_curr = state_est(3, k);
    x3_curr = state_est(5, k);
    
    v1_curr = state_est(2, k);
    v2_curr = state_est(4, k);
    v3_curr = state_est(6, k);
    
    x_vec = [x1_curr; x2_curr; x3_curr];
    v_vec = [v1_curr; v2_curr; v3_curr];
    
    % Predicted acceleration with current parameters
    a_pred = M \ (F_ext(:, k) - C_est*v_vec - K_est*x_vec);
    
    % Compute sensitivity matrix H_param (∂a/∂θ) using analytical derivatives
    H_param = zeros(3, n_params);
    
    % Sensitivity to k1 (only affects equation for mass 1)
    dK_dk1 = [1, 0, 0; 0, 0, 0; 0, 0, 0];
    H_param(:, 1) = -M \ (dK_dk1 * x_vec);
    
    % Sensitivity to k2 (affects equations for masses 1 and 2)
    dK_dk2 = [1, -1, 0; -1, 1, 0; 0, 0, 0];
    H_param(:, 2) = -M \ (dK_dk2 * x_vec);
    
    % Sensitivity to k3 (affects equations for masses 2 and 3)
    dK_dk3 = [0, 0, 0; 0, 1, -1; 0, -1, 1];
    H_param(:, 3) = -M \ (dK_dk3 * x_vec);
    
    % Sensitivity to c1
    dC_dc1 = [1, 0, 0; 0, 0, 0; 0, 0, 0];
    H_param(:, 4) = -M \ (dC_dc1 * v_vec);
    
    % Sensitivity to c2
    dC_dc2 = [1, -1, 0; -1, 1, 0; 0, 0, 0];
    H_param(:, 5) = -M \ (dC_dc2 * v_vec);
    
    % Sensitivity to c3
    dC_dc3 = [0, 0, 0; 0, 1, -1; 0, -1, 1];
    H_param(:, 6) = -M \ (dC_dc3 * v_vec);
    
    % Innovation for parameters (acceleration residuals)
    innovation_param = Z_accel - a_pred;
    
    % Measurement noise for acceleration
    R_accel = diag([noise_std_a1^2, noise_std_a2^2, noise_std_a3^2]);
    
    % Innovation covariance for parameters
    S_param = H_param * P_param_pred * H_param' + R_accel;
    
    % Regularize to avoid ill-conditioning
    S_param = S_param + 1e-6 * eye(3);
    
    % Kalman gain for parameters
    K_param = P_param_pred * H_param' / S_param;
    
    % Parameter update
    params_update = params_pred + K_param * innovation_param;
    
    % Apply physical constraints
    params_update(1:3) = max(params_update(1:3), 100);   % Minimum stiffness
    params_update(1:3) = min(params_update(1:3), 2000);  % Maximum stiffness
    params_update(4:6) = max(params_update(4:6), 0.1);   % Minimum damping
    params_update(4:6) = min(params_update(4:6), 20);    % Maximum damping
    
    % Store updated parameters
    params_est(:, k) = params_update;
    
    % Parameter covariance update
    P_param = (eye(n_params) - K_param * H_param) * P_param_pred;
    
    % Store stiffness estimates for plotting
    k1_est(k) = params_est(1, k);
    k2_est(k) = params_est(2, k);
    k3_est(k) = params_est(3, k);
end

fprintf('DEKF completed.\n');

%% ========================================================================
% VISUALIZATION - MULTIPLE DAMAGE TRACKING
% ========================================================================
fprintf('\nGenerating plots...\n');

% Figure 1: All spring stiffness tracking
figure('Position', [50, 50, 1400, 800]);

% Spring 1
subplot(3, 2, 1);
plot(time, k1_true, 'k-', 'LineWidth', 3, 'DisplayName', 'True k₁');
hold on;
plot(time, k1_est, 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated k₁');
plot(time, k1*ones(size(time)), 'b:', 'LineWidth', 1.5, 'DisplayName', 'Initial k₁');
xlabel('Time (s)', 'FontSize', 11);
ylabel('Spring Stiffness k₁ (N/m)', 'FontSize', 11);
title('Spring 1 Stiffness Tracking', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, T]);
ylim([0.5*k1, 1.1*k1]);

% Spring 2
subplot(3, 2, 3);
plot(time, k2_true, 'k-', 'LineWidth', 3, 'DisplayName', 'True k₂');
hold on;
plot(time, k2_est, 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated k₂');
plot(time, k2*ones(size(time)), 'b:', 'LineWidth', 1.5, 'DisplayName', 'Initial k₂');
xlabel('Time (s)', 'FontSize', 11);
ylabel('Spring Stiffness k₂ (N/m)', 'FontSize', 11);
title('Spring 2 Stiffness Tracking', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, T]);
ylim([0.5*k2, 1.1*k2]);

% Spring 3
subplot(3, 2, 5);
plot(time, k3_true, 'k-', 'LineWidth', 3, 'DisplayName', 'True k₃');
hold on;
plot(time, k3_est, 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated k₃');
plot(time, k3*ones(size(time)), 'b:', 'LineWidth', 1.5, 'DisplayName', 'Initial k₃');
xlabel('Time (s)', 'FontSize', 11);
ylabel('Spring Stiffness k₃ (N/m)', 'FontSize', 11);
title('Spring 3 Stiffness Tracking', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, T]);
ylim([0.5*k3, 1.1*k3]);

% Error plots
subplot(3, 2, 2);
error_k1 = abs(k1_est - k1_true) ./ k1_true * 100;
plot(time, error_k1, 'm-', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 11);
ylabel('Estimation Error (%)', 'FontSize', 11);
title('Spring 1 Estimation Error', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, T]);
ylim([0, 30]);

subplot(3, 2, 4);
error_k2 = abs(k2_est - k2_true) ./ k2_true * 100;
plot(time, error_k2, 'm-', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 11);
ylabel('Estimation Error (%)', 'FontSize', 11);
title('Spring 2 Estimation Error', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, T]);
ylim([0, 30]);

subplot(3, 2, 6);
error_k3 = abs(k3_est - k3_true) ./ k3_true * 100;
plot(time, error_k3, 'm-', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 11);
ylabel('Estimation Error (%)', 'FontSize', 11);
title('Spring 3 Estimation Error', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, T]);
ylim([0, 30]);

sgtitle('Multiple Spring Damage Tracking', 'FontSize', 16, 'FontWeight', 'bold');

% Figure 2: All stiffness tracking in one plot
figure('Position', [50, 50, 1200, 500]);

plot(time, k1_true, 'k-', 'LineWidth', 3, 'DisplayName', 'True k₁');
hold on;
plot(time, k2_true, 'k-', 'LineWidth', 3, 'DisplayName', 'True k₂');
plot(time, k3_true, 'k-', 'LineWidth', 3, 'DisplayName', 'True k₃');

plot(time, k1_est, 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated k₁');
plot(time, k2_est, 'b--', 'LineWidth', 2, 'DisplayName', 'Estimated k₂');
plot(time, k3_est, 'g--', 'LineWidth', 2, 'DisplayName', 'Estimated k₃');

% Plot initial values
plot(time, k1*ones(size(time)), 'r:', 'LineWidth', 1);
plot(time, k2*ones(size(time)), 'b:', 'LineWidth', 1);
plot(time, k3*ones(size(time)), 'g:', 'LineWidth', 1);

xlabel('Time (s)', 'FontSize', 12);
ylabel('Spring Stiffness (N/m)', 'FontSize', 12);
title('Simultaneous Damage Tracking in All Three Springs', 'FontSize', 14, 'FontWeight', 'bold');
legend({'True k₁', 'True k₂', 'True k₃', 'Estimated k₁', 'Estimated k₂', 'Estimated k₃'}, ...
       'Location', 'best', 'FontSize', 10);
grid on;
xlim([0, T]);

% Add damage percentage labels
text(1, k1*0.9, sprintf('k₁: %.0f%% damage', max(damage_k1)*100), ...
     'FontSize', 10, 'Color', 'r', 'BackgroundColor', 'white');
text(6, k2*0.9, sprintf('k₂: %.0f%% damage', max(damage_k2)*100), ...
     'FontSize', 10, 'Color', 'b', 'BackgroundColor', 'white');
text(9, k3*0.9, sprintf('k₃: %.0f%% damage', max(damage_k3)*100), ...
     'FontSize', 10, 'Color', 'g', 'BackgroundColor', 'white');

% Figure 3: Damage factor tracking (normalized)
figure('Position', [50, 50, 1200, 400]);

% True damage factors
subplot(1, 2, 1);
plot(time, damage_k1*100, 'r-', 'LineWidth', 2, 'DisplayName', 'True Damage k₁');
hold on;
plot(time, damage_k2*100, 'b-', 'LineWidth', 2, 'DisplayName', 'True Damage k₂');
plot(time, damage_k3*100, 'g-', 'LineWidth', 2, 'DisplayName', 'True Damage k₃');
xlabel('Time (s)', 'FontSize', 12);
ylabel('Damage (%)', 'FontSize', 12);
title('True Damage Progression', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, T]);
ylim([0, 70]);

% Estimated damage factors
subplot(1, 2, 2);
est_damage_k1 = 1 - k1_est/k1;
est_damage_k2 = 1 - k2_est/k2;
est_damage_k3 = 1 - k3_est/k3;

plot(time, est_damage_k1*100, 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated Damage k₁');
hold on;
plot(time, est_damage_k2*100, 'b--', 'LineWidth', 2, 'DisplayName', 'Estimated Damage k₂');
plot(time, est_damage_k3*100, 'g--', 'LineWidth', 2, 'DisplayName', 'Estimated Damage k₃');
xlabel('Time (s)', 'FontSize', 12);
ylabel('Damage (%)', 'FontSize', 12);
title('Estimated Damage Progression', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, T]);
ylim([0, 70]);

% Figure 4: State estimation performance
figure('Position', [50, 50, 1200, 600]);

% Plot position estimation for all masses
subplot(2, 1, 1);
plot(time, x1_true, 'k-', 'LineWidth', 1.5, 'DisplayName', 'True x₁');
hold on;
plot(time, x2_true, 'k-', 'LineWidth', 1.5, 'DisplayName', 'True x₂');
plot(time, x3_true, 'k-', 'LineWidth', 1.5, 'DisplayName', 'True x₃');

plot(time, state_est(1,:), 'r--', 'LineWidth', 1, 'DisplayName', 'Estimated x₁');
plot(time, state_est(3,:), 'b--', 'LineWidth', 1, 'DisplayName', 'Estimated x₂');
plot(time, state_est(5,:), 'g--', 'LineWidth', 1, 'DisplayName', 'Estimated x₃');

xlabel('Time (s)', 'FontSize', 12);
ylabel('Displacement (m)', 'FontSize', 12);
title('Position Estimation', 'FontSize', 14, 'FontWeight', 'bold');
legend({'True x₁', 'True x₂', 'True x₃', 'Estimated x₁', 'Estimated x₂', 'Estimated x₃'}, ...
       'Location', 'best', 'FontSize', 9);
grid on;
xlim([0, T]);

% Plot position errors
subplot(2, 1, 2);
error_x1 = abs(state_est(1,:) - x1_true);
error_x2 = abs(state_est(3,:) - x2_true);
error_x3 = abs(state_est(5,:) - x3_true);

plot(time, error_x1, 'r-', 'LineWidth', 1, 'DisplayName', 'Error x₁');
hold on;
plot(time, error_x2, 'b-', 'LineWidth', 1, 'DisplayName', 'Error x₂');
plot(time, error_x3, 'g-', 'LineWidth', 1, 'DisplayName', 'Error x₃');

xlabel('Time (s)', 'FontSize', 12);
ylabel('Estimation Error (m)', 'FontSize', 12);
title('Position Estimation Errors', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, T]);

% Figure 5: Performance metrics summary
figure('Position', [50, 50, 1000, 400]);

% Calculate final errors and RMSE
final_errors = zeros(3,1);
rmse_values = zeros(3,1);
for i = 1:3
    k_true = eval(sprintf('k%d_true', i));
    k_est = eval(sprintf('k%d_est', i));
    final_errors(i) = abs(k_est(end) - k_true(end)) / k_true(end) * 100;
    rmse_values(i) = sqrt(mean((k_est - k_true).^2));
end

subplot(1, 2, 1);
bar(1:3, final_errors);
set(gca, 'XTick', 1:3, 'XTickLabel', {'k₁', 'k₂', 'k₃'});
ylabel('Final Estimation Error (%)', 'FontSize', 12);
title('Final Parameter Estimation Error', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ylim([0, 20]);

% Add value labels
for i = 1:3
    text(i, final_errors(i)+0.5, sprintf('%.1f%%', final_errors(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11);
end

subplot(1, 2, 2);
bar(1:3, rmse_values);
set(gca, 'XTick', 1:3, 'XTickLabel', {'k₁', 'k₂', 'k₃'});
ylabel('RMSE (N/m)', 'FontSize', 12);
title('Root Mean Square Error', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Add value labels
for i = 1:3
    text(i, rmse_values(i)+5, sprintf('%.1f', rmse_values(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11);
end

%% ========================================================================
% PERFORMANCE ANALYSIS
% ========================================================================
fprintf('\n================================================================\n');
fprintf('PERFORMANCE ANALYSIS - MULTIPLE DAMAGE TRACKING\n');
fprintf('================================================================\n');

% RMSE calculations
k1_rmse = sqrt(mean((k1_est - k1_true).^2));
k2_rmse = sqrt(mean((k2_est - k2_true).^2));
k3_rmse = sqrt(mean((k3_est - k3_true).^2));

fprintf('\nStiffness Estimation RMSE:\n');
fprintf('  Spring 1 (k₁): %.1f N/m (%.1f%% of initial value)\n', k1_rmse, k1_rmse/k1*100);
fprintf('  Spring 2 (k₂): %.1f N/m (%.1f%% of initial value)\n', k2_rmse, k2_rmse/k2*100);
fprintf('  Spring 3 (k₃): %.1f N/m (%.1f%% of initial value)\n', k3_rmse, k3_rmse/k3*100);

% Final estimation accuracy
fprintf('\nFinal Stiffness Values (at t = %.1f s):\n', T);
fprintf('Spring\tTrue\t\tEstimated\tError\t\tDamage Detected\n');
for i = 1:3
    k_true_end = eval(sprintf('k%d_true(end)', i));
    k_est_end = eval(sprintf('k%d_est(end)', i));
    error = abs(k_est_end - k_true_end) / k_true_end * 100;
    damage_true = eval(sprintf('max(damage_k%d)*100', i));
    
    fprintf('k%d\t%.1f\t\t%.1f\t\t%.1f%%\t\t%.0f%% (True: %.0f%%)\n', ...
            i, k_true_end, k_est_end, error, (1-k_est_end/eval(sprintf('k%d', i)))*100, damage_true);
end

% State estimation performance
state_rmse = sqrt(mean((state_est - state_true).^2, 2));
fprintf('\nState Estimation RMSE:\n');
fprintf('  x₁: %.4f m\tv₁: %.4f m/s\n', state_rmse(1), state_rmse(2));
fprintf('  x₂: %.4f m\tv₂: %.4f m/s\n', state_rmse(3), state_rmse(4));
fprintf('  x₃: %.4f m\tv₃: %.4f m/s\n', state_rmse(5), state_rmse(6));

% Detection timing analysis
threshold = 0.95;  % 95% of initial value
detection_times = zeros(3,1);
true_start_times = [1, 5, 2];  % Approximate start times from damage profiles

for i = 1:3
    k_est = eval(sprintf('k%d_est', i));
    k_initial = eval(sprintf('k%d', i));
    
    % Find when stiffness drops below threshold
    idx = find(k_est < threshold * k_initial, 1);
    if ~isempty(idx)
        detection_times(i) = time(idx);
    else
        detection_times(i) = NaN;
    end
end

fprintf('\nDamage Detection Timing:\n');
fprintf('Spring\tTrue Start\tDetection\tDelay\n');
for i = 1:3
    if ~isnan(detection_times(i))
        delay = detection_times(i) - true_start_times(i);
        fprintf('k%d\t%.1f s\t\t%.1f s\t\t%.1f s\n', i, true_start_times(i), detection_times(i), delay);
    else
        fprintf('k%d\t%.1f s\t\tNot detected\t-\n', i, true_start_times(i));
    end
end

% Correlation analysis between estimated damages
fprintf('\nDamage Progression Correlation:\n');
corr_k1k2 = corr(k1_est', k2_est');
corr_k1k3 = corr(k1_est', k3_est');
corr_k2k3 = corr(k2_est', k3_est');
fprintf('  k₁ vs k₂: %.3f\n', corr_k1k2);
fprintf('  k₁ vs k₃: %.3f\n', corr_k1k3);
fprintf('  k₂ vs k₃: %.3f\n', corr_k2k3);

fprintf('\n================================================================\n');
fprintf('Simulation Complete\n');
fprintf('Total time: %.1f s, Time step: %.4f s\n', T, dt);
fprintf('Measurement noise: Pos=%.1f%%, Vel=%.1f%%, Accel=%.1f%%\n', ...
        SNR_position*100, SNR_velocity*100, SNR_acceleration*100);
fprintf('================================================================\n');
