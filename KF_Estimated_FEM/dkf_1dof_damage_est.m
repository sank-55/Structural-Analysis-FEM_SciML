% DUAL KALMAN FILTER for 1-DOF DAMAGE DETECTION
clear; clc; close all;

%% System Parameters (1 DOF)
m = 1.0;           % mass (kg)
k = 100.0;         % initial stiffness (N/m)
c = 1.5;           % damping (N·s/m)
F = 5.0;           % constant external force (N)

dt = 0.0001;         % time step (s)
T = 5;            % total time (s)
N = T / dt;        % number of time steps
t = (0:N-1) * dt;  % time vector

% True damage evolution (linear increase)
theta_true_x = linspace(0, 1, N);
theta_true = nan(1,N);
for x=1:length(theta_true_x)
    theta_true(x)=0.8*theta_true_x(x)^2*theta_true_x(x)+0.2;
end

%% Simulate True System with Damage
x_true = zeros(2, N);  % state: [displacement; velocity]
a_meas = zeros(1, N);  % measured acceleration (with noise)
noise_std = 0.0001;

for k = 1:N-1
    u = x_true(1,k); v = x_true(2,k);
    theta = theta_true(k);
    a = (1/m)*(F - c*v - (1 - theta)*k*u);  % true acceleration
    x_true(2,k+1) = v + dt * a;
    x_true(1,k+1) = u + dt * v;
    a_meas(k) = a + noise_std * randn;     % add noise
end
% last measurement
u = x_true(1,N); v = x_true(2,N);
a_meas(N) = (1/m)*(F - c*v - (1 - theta_true(N))*k*u) + noise_std * randn;

%% DUAL KALMAN FILTER INITIALIZATION
% State filter: x = [u; v]
xhat = [0; 0];               % initial state estimate
P = 0.0001*eye(2);                  % state covariance
Q = 1e-16 * eye(2);           % process noise
R = 1e-6;             % measurement noise

% Parameter filter: scalar θ
theta_hat = 0.0;             % initial damage estimate
P_theta = 0.1;               % initial variance
Q_theta = 1e-7;              % small drift allowed
R_theta = 1e-3 ;                 % assume same noise level

% Storage
theta_est = zeros(1, N);
theta_est(1) = theta_hat;

%% DUAL KF LOOP
for k = 2:N
    %-----------------------
    % 1) STATE KF
    %-----------------------
    % Use previous theta_hat to compute A
    A = [0 1;
        -(1 - theta_hat)*k/m   -c/m];
    B = [0; F/m];
    Ad = eye(2) + dt * A;
    Bd = dt * B;

    % Predict
    x_pred = Ad * xhat + Bd;
    P_pred = Ad * P * Ad' + Q;

    % Measurement: acceleration = (1/m)*(F - c*v - (1 - theta)*k*u)
    H = [- (1 - theta_hat)*k/m, -c/m];  % linear in x
    z_pred = H * x_pred + F/m;
    z_meas = a_meas(k);

    % Kalman Gain and update
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    xhat = x_pred + K * (z_meas - z_pred);
    P = (eye(2) - K * H) * P_pred;

    %-----------------------
    % 2) PARAMETER EKF
    %-----------------------
    u_hat = xhat(1);
    v_hat = xhat(2);
    % Measurement model: z = (1/m)*(F - c*v - (1 - theta)*k*u)
    % => nonlinear in theta
    z_theta = (1/m)*(F - c*v_hat - (1 - theta_hat)*k*u_hat);
    H_theta = (k/m) * u_hat;  % ∂z/∂theta

    % EKF update
    P_theta_pred = P_theta + Q_theta;
    S_theta = H_theta * P_theta_pred * H_theta' + R_theta;
    K_theta = P_theta_pred * H_theta' / S_theta;
    theta_hat = theta_hat + K_theta * (z_meas - z_theta);
    P_theta = (1 - K_theta * H_theta) * P_theta_pred;

    % Optional: constrain theta between 0 and 1
    theta_hat = max(0, min(1, theta_hat));

    theta_est(k) = theta_hat;
end

%% Plot Results
figure;
plot(t, theta_true, 'b-', 'LineWidth', 2); hold on;
plot(t, theta_est, 'r--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Damage Factor \theta');
legend('True \theta','Estimated \theta');
title('1-DOF Damage Detection Using Dual Kalman Filter');
grid on;
