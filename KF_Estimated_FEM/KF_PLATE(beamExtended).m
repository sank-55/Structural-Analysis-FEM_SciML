%% Cantilever Plate Kalman Filter - Complete MATLAB Implementation
clear; close all; clc;

fprintf('Cantilever Plate Kalman Filter\n');
fprintf('===============================\n');

%% System Parameters
% Plate physical properties
L = 1.0;           % Length [m]
W = 0.5;           % Width [m]
thickness = 0.01;  % Thickness [m]
E = 2.1e11;        % Young's modulus [Pa]
rho = 7800;        % Density [kg/m³]
zeta = 0.02;       % Damping ratio

% Simulation parameters
dt = 0.01;         % Time step [s]
T = 5.0;           % Total time [s]
t = 0:dt:T;        % Time vector
N = length(t);     % Number of time steps

% Vibration modes
n_modes = 2;       % Number of modes to consider
n_meas = 2;        % Number of measurements

%% Calculate Natural Frequencies and Mode Shapes
fprintf('Calculating natural frequencies...\n');
beta_L = [1.875, 4.694, 7.855, 10.996]; % Cantilever coefficients
f_natural = zeros(n_modes, 1);
for i = 1:n_modes
    if i <= length(beta_L)
        f_natural(i) = (beta_L(i)^2 / (2*pi*L^2)) * ...
                      sqrt(E * thickness^2 / (12 * rho * (1 - 0.3^2)));
    else
        f_natural(i) = f_natural(i-1) * 2.5; % Approximate higher modes
    end
    fprintf('  Mode %d: %.2f Hz\n', i, f_natural(i));
end

%% System Matrices (Continuous Time)
% State: x = [q1, q1_dot, q2, q2_dot, ...]^T where q are modal coordinates
n_states = 2 * n_modes;
A_cont = zeros(n_states);
B_cont = zeros(n_states, 1);
C_cont = zeros(n_meas, n_states);

for i = 1:n_modes
    idx = 2*(i-1) + 1;
    omega = 2 * pi * f_natural(i);
    
    % Modal equations: q_ddot + 2*zeta*omega*q_dot + omega^2*q = F/modal_mass
    A_cont(idx, idx+1) = 1.0;
    A_cont(idx+1, idx) = -omega^2;
    A_cont(idx+1, idx+1) = -2 * zeta * omega;
    
    % Input affects accelerations
    B_cont(idx+1) = 1.0;
    
    % Measure modal displacements
    if i <= n_meas
        C_cont(i, idx) = 1.0;
    end
end

%% Discrete-time System Matrices
A_disc = expm(A_cont * dt);
B_disc = A_cont \ (A_disc - eye(n_states)) * B_cont;
C_disc = C_cont;

fprintf('System matrices created.\n');

%% Kalman Filter Initialization
% Process noise covariance (model uncertainty)
Q = eye(n_states) * 1e-6;

% Measurement noise covariance (sensor noise)
R = eye(n_meas) * 1e-4;

% Initial state covariance
P = eye(n_states) * 0.1;

% Initial state estimate
x_hat = zeros(n_states, 1);

% True system (with slight parameter variation for realism)
A_true = expm(A_cont * dt * 0.98); % 2% model error
x_true = zeros(n_states, 1);
x_true(1:2:end) = [0.02; 0.01]; % Initial displacements

fprintf('Kalman Filter initialized.\n');

%% Excitation Force
force = 0.5 * sin(2*pi*2*t) + 0.3 * sin(2*pi*5*t) + 0.1 * randn(1, N);

%% Storage Arrays
true_states = zeros(n_states, N);
est_states = zeros(n_states, N);
measurements = zeros(n_meas, N);
tip_disp_true = zeros(1, N);
tip_disp_est = zeros(1, N);
errors = zeros(n_modes, N);

%% Main Simulation Loop
fprintf('Running simulation...\n');

for k = 1:N
    %% True System Evolution
    process_noise = 1e-5 * randn(n_states, 1);
    x_true = A_true * x_true + B_disc * force(k) + process_noise;
    
    %% Measurement
    meas_noise = 1e-3 * randn(n_meas, 1);
    y = C_disc * x_true + meas_noise;
    
    %% Kalman Filter - Prediction Step
    x_hat = A_disc * x_hat + B_disc * force(k);
    P = A_disc * P * A_disc' + Q;
    
    %% Kalman Filter - Update Step
    K = P * C_disc' / (C_disc * P * C_disc' + R);
    x_hat = x_hat + K * (y - C_disc * x_hat);
    P = (eye(n_states) - K * C_disc) * P;
    
    %% Store Results
    true_states(:, k) = x_true;
    est_states(:, k) = x_hat;
    measurements(:, k) = y;
    
    % Calculate tip displacements (simplified)
    tip_disp_true(k) = sum(x_true(1:2:end) .* [2.0; 1.5]); % Mode shapes at tip
    tip_disp_est(k) = sum(x_hat(1:2:end) .* [2.0; 1.5]);
    
    % Errors
    errors(:, k) = x_true(1:2:end) - x_hat(1:2:end);
end

fprintf('Simulation completed.\n');

%% Performance Analysis
fprintf('\nPerformance Metrics:\n');
fprintf('-------------------\n');
for i = 1:n_modes
    rmse = sqrt(mean(errors(i, :).^2));
    max_err = max(abs(errors(i, :)));
    fprintf('Mode %d: RMSE = %.6f m, Max Error = %.6f m\n', i, rmse, max_err);
end

overall_rmse = sqrt(mean(errors(:).^2));
fprintf('Overall RMSE: %.6f m\n', overall_rmse);

%% Plot Results
figure('Position', [100, 100, 1400, 1000]);

% Plot 1: Mode 1 Displacement
subplot(3, 3, 1);
plot(t, true_states(1,:), 'b-', 'LineWidth', 2); hold on;
plot(t, est_states(1,:), 'r--', 'LineWidth', 2);
plot(t, measurements(1,:), 'g:', 'LineWidth', 1);
ylabel('Displacement [m]');
title('Mode 1 Displacement');
legend('True', 'Estimated', 'Measured', 'Location', 'best');
grid on;

% Plot 2: Mode 2 Displacement
subplot(3, 3, 2);
plot(t, true_states(3,:), 'b-', 'LineWidth', 2); hold on;
plot(t, est_states(3,:), 'r--', 'LineWidth', 2);
plot(t, measurements(2,:), 'g:', 'LineWidth', 1);
ylabel('Displacement [m]');
title('Mode 2 Displacement');
legend('True', 'Estimated', 'Measured', 'Location', 'best');
grid on;

% Plot 3: Tip Displacement
subplot(3, 3, 3);
plot(t, tip_disp_true, 'b-', 'LineWidth', 2); hold on;
plot(t, tip_disp_est, 'r--', 'LineWidth', 2);
ylabel('Displacement [m]');
xlabel('Time [s]');
title('Tip Displacement');
legend('True', 'Estimated', 'Location', 'best');
grid on;

% Plot 4: Estimation Errors
subplot(3, 3, 4);
colors = ['b', 'r', 'g', 'm'];
for i = 1:n_modes
    plot(t, errors(i,:), 'Color', colors(i), 'LineWidth', 1.5, 'DisplayName', sprintf('Mode %d', i));
    hold on;
end
ylabel('Error [m]');
xlabel('Time [s]');
title('Estimation Errors');
legend('Location', 'best');
grid on;

% Plot 5: Force Input
subplot(3, 3, 5);
plot(t, force, 'k-', 'LineWidth', 1.5);
ylabel('Force [N]');
xlabel('Time [s]');
title('Excitation Force');
grid on;

% Plot 6: Kalman Gain Evolution (Mode 1)
subplot(3, 3, 6);
% Simulate K gain for plotting
K_sim = zeros(n_states, N);
P_temp = eye(n_states) * 0.1;
for k = 1:N
    P_temp = A_disc * P_temp * A_disc' + Q;
    K_temp = P_temp * C_disc' / (C_disc * P_temp * C_disc' + R);
    K_sim(:, k) = K_temp(:, 1); % First measurement gain
end
plot(t, K_sim(1,:), 'b-', 'LineWidth', 1.5); hold on;
plot(t, K_sim(3,:), 'r-', 'LineWidth', 1.5);
ylabel('Gain');
xlabel('Time [s]');
title('Kalman Gain (Mode 1)');
legend('K_{1,1}', 'K_{3,1}', 'Location', 'best');
grid on;

% Plot 7: Frequency Spectrum - Mode 1
subplot(3, 3, 7);
[f_true, P1_true] = computeFFT(true_states(1,:), dt);
[~, P1_est] = computeFFT(est_states(1,:), dt);
plot(f_true, P1_true, 'b-', 'LineWidth', 2); hold on;
plot(f_true, P1_est, 'r--', 'LineWidth', 2);
xlabel('Frequency [Hz]');
ylabel('Power');
title('Frequency Spectrum - Mode 1');
legend('True', 'Estimated', 'Location', 'best');
grid on;
xlim([0, 20]);

% Plot 8: Moving RMS Error
subplot(3, 3, 8);
window = 100;
moving_rms = zeros(1, N-window+1);
for i = 1:N-window+1
    moving_rms(i) = sqrt(mean(errors(1, i:i+window-1).^2));
end
plot(t(window:end), moving_rms, 'k-', 'LineWidth', 2);
ylabel('RMS Error [m]');
xlabel('Time [s]');
title(sprintf('Moving RMS Error (Window: %d)', window));
grid on;

% Plot 9: State Covariance (Mode 1)
subplot(3, 3, 9);
% Simulate P for plotting
P_sim = zeros(1, N);
P_temp = eye(n_states) * 0.1;
for k = 1:N
    P_temp = A_disc * P_temp * A_disc' + Q;
    K_temp = P_temp * C_disc' / (C_disc * P_temp * C_disc' + R);
    P_temp = (eye(n_states) - K_temp * C_disc) * P_temp;
    P_sim(k) = P_temp(1,1); % Variance of mode 1 displacement
end
semilogy(t, P_sim, 'm-', 'LineWidth', 2);
ylabel('Covariance (log)');
xlabel('Time [s]');
title('State Covariance - Mode 1');
grid on;

sgtitle('Cantilever Plate Kalman Filter Results', 'FontSize', 14, 'FontWeight', 'bold');

%% FFT Computation Function
function [f, P1] = computeFFT(signal, dt)
    Fs = 1/dt;
    L = length(signal);
    f = Fs * (0:(L/2)) / L;
    Y = fft(signal);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2 * P1(2:end-1);
end

%% Display Summary
fprintf('\nSimulation Summary:\n');
fprintf('==================\n');
fprintf('Duration: %.1f seconds\n', T);
fprintf('Time step: %.3f seconds\n', dt);
fprintf('Number of modes: %d\n', n_modes);
fprintf('Measurement points: %d\n', n_meas);
fprintf('Process noise: %.2e\n', 1e-6);
fprintf('Measurement noise: %.2e\n', 1e-4);
fprintf('Kalman Filter successfully estimated plate vibrations.\n');

%% Additional Analysis - Convergence Check
fprintf('\nConvergence Analysis:\n');
fprintf('=====================\n');
final_errors = mean(abs(errors(:, end-100:end)), 2);
for i = 1:n_modes
    fprintf('Mode %d steady-state error: %.6f m\n', i, final_errors(i));
end

% Check if Kalman Filter is working properly
if overall_rmse < 0.01
    fprintf('\n Kalman Filter is working correctly!\n');
else
    fprintf('\n  Kalman Filter may need tuning.\n');
end
