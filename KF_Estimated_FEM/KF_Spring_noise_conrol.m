%%  Kalman filter in 1D spring operation (with controllable measurement sparsity)
clear; close all; clc;
rng(0);

%% ---------------- System parameters (continuous) ----------------
m = 1.0;        % kg
k = 20.0;       % N/m
c = 0.5;        % N*s/m

% State U = [u; v]
A = [0 1; -k/m -c/m];
B = [0; 1/m];
H = [1 0];      % measure displacement only

%% ---------------- Time discretization for KF + simulation ----------------
dt = 0.01;         % discrete time step for KF and simulated measurements
T = 20;            % total time (s)
t_disc = 0:dt:T;   % discrete time grid
N = length(t_disc);

%% ---------------- Forcing (same for theory & KF simulation) -------------
f_amp = 1.0; f_freq = 0.5;  % Hz
zfun = @(tt) f_amp * sin(2*pi*f_freq * tt);

%% ---------------- Theoretical continuous solution (ODE) ----------------
odefun = @(t,x) [ x(2);
                 (1/m)*( zfun(t) - c*x(2) - k*x(1) ) ];
u0 = 0; v0 = 0;
x0 = [u0; v0];
t_fine = linspace(0, T, 5000);
[tt_fine, x_fine] = ode45(odefun, t_fine, x0);
u_theory = x_fine(:,1);
u_theory_disc = interp1(tt_fine, u_theory, t_disc);

%% ---------------- Discretize (exact ZOH) for KF ----------------
Maug = [A, B; zeros(1,3)];
expM = expm(Maug*dt);
F = expM(1:2,1:2);
G = expM(1:2,3);

%% ---------------- Process & measurement noise ----------------
Q = diag([1e-6, 1e-4]);     % process noise covariance (tune)
sigma_y = 0.3;             % measurement noise std (m)
R = sigma_y^2;

Q = (Q + Q')/2;
jitter = 1e-12;
LQ = chol(Q + jitter*eye(size(Q)), 'lower');

%% ---------------- Simulate "true" discrete system (no measurements yet) ----------------
U_true_disc = zeros(2, N);
U_true_disc(:,1) = [0; 0];
for k = 1:N-1
    wk = LQ * randn(2,1);
    U_true_disc(:,k+1) = F * U_true_disc(:,k) + G * zfun(t_disc(k)) + wk;
end

%% ---------------- Measurement sparsity control (choose one) ----------------
% Options:
%   measurement_mode = 'regular' -> measure every downsample_factor steps
%   measurement_mode = 'random'  -> measure a random subset with fraction measurement_fraction
measurement_mode = 'regular';   % 'regular' or 'random'
downsample_factor   = 200;        % used when measurement_mode == 'regular'  // more value less measurment pt,  vice verca
measurement_fraction = 0.20;    % used when measurement_mode == 'random' (0..1)
random_seed = 1;                % seed for reproducibility of random selection

switch lower(measurement_mode)
    case 'regular'
        measurement_indices = 1:downsample_factor:N;
    case 'random'
        rng(random_seed);
        num_meas = max(1, round(measurement_fraction * N));
        measurement_indices = sort(randperm(N, num_meas));
    otherwise
        error('Unknown measurement_mode. Use ''regular'' or ''random''.');
end

% Build measurement vector y: NaN where no measurement
y = nan(1, N);
for idx = measurement_indices
    y(idx) = H * U_true_disc(:, idx) + sigma_y * randn;
end

%% ---------------- Kalman Filter (predict-only at no-measurement steps) ----------------
U_est = zeros(2, N);
P_store = zeros(2,2,N);
U_est(:,1) = [0.5; 0.0];                 % deliberately off to show convergence
P_store(:,:,1) = diag([0.5, 0.5]);
K_store = zeros(2, N);
innov = nan(1, N);
S_store = nan(1, N);

for k = 1:N-1
    % Prediction
    Ubar = F * U_est(:,k) + G * zfun(t_disc(k));
    Pbar = F * P_store(:,:,k) * F' + Q;
    
    % Check whether we have a measurement at time k+1
    if ismember(k+1, measurement_indices)
        S = H * Pbar * H' + R;
        K = (Pbar * H') / S;
        yk = y(k+1);
        innov_k = yk - H * Ubar;
        Upost = Ubar + K * innov_k;
        Ppost = Pbar - K * H * Pbar;
    else
        % no measurement --> predict-only (no update)
        K = zeros(2,1);
        S = NaN;
        innov_k = NaN;
        Upost = Ubar;
        Ppost = Pbar;
    end
    
    % store
    U_est(:,k+1) = Upost;
    P_store(:,:,k+1) = Ppost;
    K_store(:,k+1) = K;
    innov(k+1) = innov_k;
    S_store(k+1) = S;
end

%% ---------------- Compute simple RMSEs ----------------
pos_rmse_vs_theory = sqrt(mean((U_est(1,:) - u_theory_disc).^2));
vel_rmse_vs_theory = sqrt(mean((U_est(2,:) - interp1(tt_fine, x_fine(:,2), t_disc)).^2));
fprintf('RMSE (KF estimate vs continuous theory): pos = %.4e m, vel = %.4e m/s\n', ...
        pos_rmse_vs_theory, vel_rmse_vs_theory);

% RMSE of measurements vs theory (only at measured times)
meas_idx = measurement_indices;
pos_rmse_meas = sqrt(mean((y(meas_idx) - u_theory_disc(meas_idx)).^2));
fprintf('RMSE (measurements vs theory) using %d measurements: pos = %.4e m\n', ...
        numel(meas_idx), pos_rmse_meas);

%% ---------------- Main focused plot: theory (ODE) vs measurements vs KF -----------
figure('Name','Theory (ODE) vs KF & measurements','Position',[100 100 900 420]);

% theoretical continuous solution (smooth)
h1 = plot(tt_fine, u_theory, 'k-', 'LineWidth', 1.6); hold on;

% measurements (only at measurement indices)
h2 = plot(t_disc(meas_idx), y(meas_idx), '.', 'Color', [0.8 0 0], 'MarkerSize', 8);

% KF posterior estimate (discrete samples, dashed)
h3 = plot(t_disc, U_est(1,:), 'b--', 'LineWidth', 1.4);

% ±2σ band from posterior P_store
pos_std = squeeze(sqrt(squeeze(P_store(1,1,:))))'; % 1 x N
u_est_fine = interp1(t_disc, U_est(1,:), tt_fine, 'linear');
pos_std_fine = interp1(t_disc, pos_std, tt_fine, 'linear');
upper_fine = u_est_fine + 2*pos_std_fine;
lower_fine = u_est_fine - 2*pos_std_fine;
hfill = fill([tt_fine fliplr(tt_fine)], [upper_fine fliplr(lower_fine)], 'b');
set(hfill, 'FaceAlpha', 0.12, 'EdgeColor', 'none');

xlabel('Time (s)'); ylabel('Displacement u (m)');
title(sprintf('Theory vs noisy measurements (N_{meas}=%d) vs KF estimate', numel(meas_idx)));
legend([h1 h2 h3], {'Theory (ODE45)','Measurements','KF estimate (u_{post})'}, 'Location','best');
grid on; xlim([0 T]);

% annotate RMSEs & measurement info
txt = sprintf('RMSE_{KF vs theory} = %.3e m\nRMSE_{meas vs theory} = %.3e m\nMeasurements = %d of %d', ...
              pos_rmse_vs_theory, pos_rmse_meas, numel(meas_idx), N);
xloc = 0.02*T; yloc = min(u_theory) + 0.8*(max(u_theory)-min(u_theory));
text(xloc, yloc, txt, 'BackgroundColor','w', 'EdgeColor','k');

%% ---------------- Compact diagnostics: innovation plot & Kalman gain -----------
figure('Name','KF diagnostics compact','Position',[150 200 700 300]);
subplot(2,1,1);
plot(t_disc, K_store(1,:), 'LineWidth', 1.2); hold on;
plot(t_disc, K_store(2,:), '--', 'LineWidth', 1.2);
legend('K(1)','K(2)'); ylabel('Kalman gain'); grid on; title('Kalman gain components');

subplot(2,1,2);
plot(t_disc, innov, 'k-', 'LineWidth', 1);
hold on;
plot(t_disc, sqrt(S_store), 'r--', t_disc, -sqrt(S_store), 'r--');
xlabel('Time (s)'); ylabel('Innovation'); title('Innovation and ±1σ'); grid on;

%% ---------------- Minimal animation (optional) ----------------
figure('Name','Small animation','Position',[400 100 700 140]);
axis equal; xlim([-1.5 1.5]); ylim([-0.5 0.5]); hold on;
plot([-1.5 1.5],[0 0],'k-','LineWidth',2);
h_true = plot(U_true_disc(1,1), 0, 'ko', 'MarkerSize', 14, 'MarkerFaceColor','k');
h_kf   = plot(U_est(1,1), 0.12, 'bo', 'MarkerSize', 10, 'MarkerFaceColor','b');
legend('rail','True discrete','KF estimate','Location','northwest');
title('True discrete vs KF estimate');

stepSkip = max(1, round(N/300));
for k=1:stepSkip:N
    set(h_true,'XData', U_true_disc(1,k));
    set(h_kf,  'XData', U_est(1,k));
    drawnow;
end

fprintf('Plotting complete. The black line is the continuous ODE solution (no process noise).\n');

    