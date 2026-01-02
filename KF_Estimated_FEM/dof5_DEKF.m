%% ----------------------------------------------------------------------------------------------------
% This code properly demonstate the Use of Dual KF system for damage detection in 5 dof system 
%% ======================================================

clear; close all; clc;

%% ========================================================================
% SYSTEM DEFINITION - 5-DOF MASS-SPRING-DAMPER
% ========================================================================
fprintf('=== 5-DOF Mass-Spring-Damper System with Multiple Damages ===\n');

% True system parameters
m1 = 10; m2 = 8; m3 = 6; m4 = 5; m5 = 4;       % Masses (kg)
k1 = 1000; k2 = 800; k3 = 600; k4 = 500; k5 = 400; % Initial spring stiffness (N/m)
c1 = 8; c2 = 6; c3 = 4; c4 = 3; c5 = 2;        % Damping coefficients (Ns/m)

% Mass matrix
M = diag([m1, m2, m3, m4, m5]);

% Initial stiffness matrix (5x5 tridiagonal)
K_initial = zeros(5,5);
K_initial(1,1) = k1 + k2;
K_initial(1,2) = -k2;
K_initial(2,1) = -k2;
K_initial(2,2) = k2 + k3;
K_initial(2,3) = -k3;
K_initial(3,2) = -k3;
K_initial(3,3) = k3 + k4;
K_initial(3,4) = -k4;
K_initial(4,3) = -k4;
K_initial(4,4) = k4 + k5;
K_initial(4,5) = -k5;
K_initial(5,4) = -k5;
K_initial(5,5) = k5;

% Initial damping matrix (proportional damping, tridiagonal)
C_initial = zeros(5,5);
C_initial(1,1) = c1 + c2;
C_initial(1,2) = -c2;
C_initial(2,1) = -c2;
C_initial(2,2) = c2 + c3;
C_initial(2,3) = -c3;
C_initial(3,2) = -c3;
C_initial(3,3) = c3 + c4;
C_initial(3,4) = -c4;
C_initial(4,3) = -c4;
C_initial(4,4) = c4 + c5;
C_initial(4,5) = -c5;
C_initial(5,4) = -c5;
C_initial(5,5) = c5;

fprintf('Initial parameters:\n');
fprintf('k1 = %.0f N/m, k2 = %.0f N/m, k3 = %.0f N/m, k4 = %.0f N/m, k5 = %.0f N/m\n', ...
        k1, k2, k3, k4, k5);
fprintf('c1 = %.1f Ns/m, c2 = %.1f Ns/m, c3 = %.1f Ns/m, c4 = %.1f Ns/m, c5 = %.1f Ns/m\n', ...
        c1, c2, c3, c4, c5);

%% ========================================================================
% TRUE DAMAGE SCENARIO - ALL SPRINGS DAMAGED PROGRESSIVELY
% ========================================================================
% Time parameters
dt = 0.005;          % Time step
T = 20;              % Total simulation time (s)
N = round(T/dt);     % Number of time steps
time = linspace(0, T, N+1);

%% Generating True Damage for 5 springs
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

% Spring 4: Late damage
damage_k4 = zeros(1, N+1);
for i = 1:N+1
    if time(i) <= 6
        damage_k4(i) = 0;  % No damage initially
    elseif time(i) <= 15
        damage_k4(i) = 0.35 * ((time(i)-6)/9);  % 35% damage by 15s
    else
        damage_k4(i) = 0.35;  % Total 35% damage
    end
end

% Spring 5: Intermittent damage
damage_k5 = zeros(1, N+1);
for i = 1:N+1
    if time(i) <= 4
        damage_k5(i) = 0.2 * (time(i)/4);  % 20% damage by 4s
    elseif time(i) <= 10
        damage_k5(i) = 0.2;  % Constant
    elseif time(i) <= 16
        damage_k5(i) = 0.2 + 0.25 * ((time(i)-10)/6); % Additional 25% by 16s
    else
        damage_k5(i) = 0.45;  % Total 45% damage
    end
end

% True stiffness values over time
k1_true = k1 * (1 - damage_k1);
k2_true = k2 * (1 - damage_k2);
k3_true = k3 * (1 - damage_k3);
k4_true = k4 * (1 - damage_k4);
k5_true = k5 * (1 - damage_k5);

fprintf('\nDamage scenario:\n');
fprintf('Spring 1: Progressive damage up to %.0f%%\n', max(damage_k1)*100);
fprintf('Spring 2: Progressive damage up to %.0f%%\n', max(damage_k2)*100);
fprintf('Spring 3: Progressive damage up to %.0f%%\n', max(damage_k3)*100);
fprintf('Spring 4: Progressive damage up to %.0f%%\n', max(damage_k4)*100);
fprintf('Spring 5: Progressive damage up to %.0f%%\n', max(damage_k5)*100);

%% ========================================================================
% SIMULATE TRUE SYSTEM RESPONSE WITH PROGRESSIVE DAMAGE
% ========================================================================
fprintf('\nSimulating true system response with progressive damage...\n');

% State vector: [x1; v1; x2; v2; x3; v3; x4; v4; x5; v5]
state_true = zeros(10, N+1);

% External force
rng(42);  % For reproducibility
F_ext = zeros(5, N+1);
for i = 1:N+1
    F_ext(:,i) = 0.1; % External force of 2N
end

% Time integration using Newmark-beta
beta_nm = 0.25;
gamma_nm = 0.5;

% Initialize state
state_true(:,1) = [0.00; 0; 0.00; 0; 0.000; 0; 0.000; 0; 0.000; 0];

for i = 1:N
    % Current stiffness values
    k1_curr = k1_true(i);
    k2_curr = k2_true(i);
    k3_curr = k3_true(i);
    k4_curr = k4_true(i);
    k5_curr = k5_true(i);
    
    % Current stiffness matrix
    K_curr = zeros(5,5);
    K_curr(1,1) = k1_curr + k2_curr;
    K_curr(1,2) = -k2_curr;
    K_curr(2,1) = -k2_curr;
    K_curr(2,2) = k2_curr + k3_curr;
    K_curr(2,3) = -k3_curr;
    K_curr(3,2) = -k3_curr;
    K_curr(3,3) = k3_curr + k4_curr;
    K_curr(3,4) = -k4_curr;
    K_curr(4,3) = -k4_curr;
    K_curr(4,4) = k4_curr + k5_curr;
    K_curr(4,5) = -k5_curr;
    K_curr(5,4) = -k5_curr;
    K_curr(5,5) = k5_curr;
    
    % Current state
    x_curr = state_true(1:2:10, i);
    v_curr = state_true(2:2:10, i);
    
    % Compute acceleration
    a_curr = M \ (F_ext(:,i) - C_initial*v_curr - K_curr*x_curr);
    
    % Predict displacement and velocity
    x_pred = x_curr + dt*v_curr + (0.5-beta_nm)*dt^2*a_curr;
    v_pred = v_curr + (1-gamma_nm)*dt*a_curr;
    
    % Next stiffness values
    k1_next = k1_true(i+1);
    k2_next = k2_true(i+1);
    k3_next = k3_true(i+1);
    k4_next = k4_true(i+1);
    k5_next = k5_true(i+1);
    
    K_next = zeros(5,5);
    K_next(1,1) = k1_next + k2_next;
    K_next(1,2) = -k2_next;
    K_next(2,1) = -k2_next;
    K_next(2,2) = k2_next + k3_next;
    K_next(2,3) = -k3_next;
    K_next(3,2) = -k3_next;
    K_next(3,3) = k3_next + k4_next;
    K_next(3,4) = -k4_next;
    K_next(4,3) = -k4_next;
    K_next(4,4) = k4_next + k5_next;
    K_next(4,5) = -k5_next;
    K_next(5,4) = -k5_next;
    K_next(5,5) = k5_next;
    
    % Solve for next acceleration
    a_next = M \ (F_ext(:,i+1) - C_initial*v_pred - K_next*x_pred);
    
    % Update velocity and displacement
    v_next = v_pred + gamma_nm*dt*a_next;
    x_next = x_pred + beta_nm*dt^2*a_next;
    
    % Store
    state_true(:, i+1) = [x_next(1); v_next(1); x_next(2); v_next(2); ...
                         x_next(3); v_next(3); x_next(4); v_next(4); ...
                         x_next(5); v_next(5)];
end

%% ========================================================================
% GENERATE NOISY MEASUREMENTS
% ========================================================================
fprintf('Generating noisy measurements...\n');

% Measurement noise levels
noise_position = 0.002;   % 0.2% noise for positions
noise_velocity = 0.003;   % 0.3% noise for velocities
noise_acceleration = 0.004; % 0.4% noise for accelerations

% Extract true states
x_true = state_true(1:2:10, :);  % Positions
v_true = state_true(2:2:10, :);  % Velocities

% Calculate true accelerations
a_true = zeros(5, N+1);
for i = 1:N+1
    % Current stiffness values
    k1_curr = k1_true(i);
    k2_curr = k2_true(i);
    k3_curr = k3_true(i);
    k4_curr = k4_true(i);
    k5_curr = k5_true(i);
    
    % Current stiffness matrix
    K_curr = zeros(5,5);
    K_curr(1,1) = k1_curr + k2_curr;
    K_curr(1,2) = -k2_curr;
    K_curr(2,1) = -k2_curr;
    K_curr(2,2) = k2_curr + k3_curr;
    K_curr(2,3) = -k3_curr;
    K_curr(3,2) = -k3_curr;
    K_curr(3,3) = k3_curr + k4_curr;
    K_curr(3,4) = -k4_curr;
    K_curr(4,3) = -k4_curr;
    K_curr(4,4) = k4_curr + k5_curr;
    K_curr(4,5) = -k5_curr;
    K_curr(5,4) = -k5_curr;
    K_curr(5,5) = k5_curr;
    
    x_curr = x_true(:, i);
    v_curr = v_true(:, i);
    a_true(:, i) = M \ (F_ext(:,i) - C_initial*v_curr - K_curr*x_curr);
end

% Add noise to measurements
x_meas = zeros(5, N+1);
v_meas = zeros(5, N+1);
a_meas = zeros(5, N+1);

for i = 1:5
    noise_std_x = sqrt(mean(x_true(i,:).^2)) * noise_position;
    noise_std_v = sqrt(mean(v_true(i,:).^2)) * noise_velocity;
    noise_std_a = sqrt(mean(a_true(i,:).^2)) * noise_acceleration;
    
    x_meas(i,:) = x_true(i,:) + noise_std_x * randn(size(x_true(i,:)));
    v_meas(i,:) = v_true(i,:) + noise_std_v * randn(size(v_true(i,:)));
    a_meas(i,:) = a_true(i,:) + noise_std_a * randn(size(a_true(i,:)));
end

fprintf('Measurement noise added:\n');
fprintf('  Positions: %.1f%%\n', noise_position*100);
fprintf('  Velocities: %.1f%%\n', noise_velocity*100);
fprintf('  Accelerations: %.1f%%\n', noise_acceleration*100);

%% ========================================================================
% DUAL EXTENDED KALMAN FILTER - MULTIPLE DAMAGE TRACKING (5 DOF)
% ========================================================================
fprintf('\n=== DUAL EXTENDED KALMAN FILTER - MULTIPLE DAMAGE TRACKING ===\n');

% Parameters to estimate: [k1, k2, k3, k4, k5, c1, c2, c3, c4, c5]
n_params = 10;  % 5 stiffness + 5 damping
n_states = 10;  % [x1, v1, x2, v2, x3, v3, x4, v4, x5, v5]

% Initial parameter estimates (slightly different from true)
params_est = zeros(n_params, N+1);
params_est(:,1) = [k1*1; k2*1; k3*1.0; k4*0.99; k5*1.01; ...
                   c1*0.998; c2*1.002; c3*0.999; c4*1.001; c5*1];

% Initial damage estimation
damage_est = zeros(5, N+1);  % Only stiffness damage

% State estimates (two approaches)
state_est = zeros(n_states, N+1);  % ki as parameters
state_est_damage = zeros(n_states, N+1);  % damage factor as parameters

% Initialize states
state_est(:,1) = state_true(:,1);
state_est_damage(:,1) = state_true(:,1);

% Covariance matrices
P_state = eye(n_states) * 1e-2;
P_param = eye(n_params) * 1e-1;
P_state_damage = P_param;
P_param_damage = 1e-1 *eye(5); % as fd is 5 for K only

% Process noise covariance
Q_state = diag([1e-4, 1e-6, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4]);
Q_param = diag([1e-1*ones(1,5), 1e-1*ones(1,5)]);  % More noise for stiffness
Q_param_damage = diag(1e-1 *ones(1,5));

% Measurement noise covariance
R_state = diag([noise_position^2 * mean(x_true(1,:).^2), ...
                noise_velocity^2 * mean(v_true(1,:).^2), ...
                noise_position^2 * mean(x_true(2,:).^2), ...
                noise_velocity^2 * mean(v_true(2,:).^2), ...
                noise_position^2 * mean(x_true(3,:).^2), ...
                noise_velocity^2 * mean(v_true(3,:).^2), ...
                noise_position^2 * mean(x_true(4,:).^2), ...
                noise_velocity^2 * mean(v_true(4,:).^2), ...
                noise_position^2 * mean(x_true(5,:).^2), ...
                noise_velocity^2 * mean(v_true(5,:).^2)]);

% Measurement matrix
H_state = eye(n_states);

% Store stiffness estimates
k_est = zeros(5, N+1);  % ki as parameters
k_est_damage = zeros(5, N+1);  % damage factor as parameters
k_est_damage_curr = zeros(5,1); % storing the k value at that moment of time 
% to show in viz
k_est_damage_ = zeros(5,N+1);
k_est_damage_(:,1)=[k1;k2;k3;k4;k5]; 

for i = 1:5
    k_est(i,1) = params_est(i,1);
    k_est_damage(i,1) = eval(sprintf('k%d', i)) * (1 - damage_est(i,1));
end

% Progress indicator
progress_step = round(N/20);
fprintf('Running DEKF for multiple damage tracking...\n');

% Forgetting factor
lambda = 1;

for k = 2:N+1
    if mod(k, progress_step) == 0
        fprintf('  Progress: %d%%\n', round(100*k/(N+1)));
    end
    
    % ====================================================================
    % APPROACH 1: ki as parameters
    % ====================================================================
    
    % Current parameter estimates
    k_est_curr = params_est(1:5, k-1);
    c_est = params_est(6:10, k-1);

    % Building Kest for damage factor;
    for i = 1:5
        k_est_damage_curr(i) = eval(sprintf('k%d', i)) * (1 - damage_est(i, k-1));
    end
    k_est_damage=build_stiffness_matrix(k_est_damage_curr);
    
    % Build stiffness and damping matrices
    K_est = build_stiffness_matrix(k_est_curr);
    C_est = build_damping_matrix(c_est);
    
    % State prediction
    A_cont = [zeros(5), eye(5); -M\K_est, -M\C_est];
    B_cont = [zeros(5); inv(M)];
    A_disc = expm(A_cont * dt);
    B_disc = (A_cont\(A_disc - eye(10))) * B_cont;

    state_pred = A_disc * state_est(:, k-1) + B_disc * F_ext(:, k-1);
    P_state_pred = A_disc * P_state * A_disc' / lambda + Q_state;

    % For damage factor 
    A_cont_damage = [zeros(5), eye(5); -M\k_est_damage, -M\C_est];
    B_cont_damage = [zeros(5); inv(M)];
    A_disc_damage = expm(A_cont_damage * dt);
    B_disc_damage = (A_cont_damage\(A_disc_damage - eye(10))) * B_cont_damage;
    
    state_pred_damage = A_disc_damage * state_est_damage(:, k-1) + B_disc_damage * F_ext(:, k-1);
    P_state_pred_damage = A_disc_damage * P_state_damage * A_disc_damage' / lambda + Q_state;
    
 
    Z_state = reshape([x_meas(:,k)'; v_meas(:,k)'], [], 1);
     % State update for ki as param
    innovation_state = Z_state - H_state * state_pred;
    S_state = H_state * P_state_pred * H_state' + R_state;
    K_state = P_state_pred * H_state' / S_state;
    state_est(:, k) = state_pred + K_state * innovation_state;
    P_state = (eye(n_states) - K_state * H_state) * P_state_pred;
    
    % State update for ki as param
    innovation_state_ = Z_state - H_state * state_pred_damage;
    S_state_ = H_state * P_state_pred_damage * H_state' + R_state;
    K_state_ = P_state_pred_damage * H_state' / S_state_;
    state_est_damage(:, k) = state_pred_damage + K_state_ * innovation_state_;
    P_state_damage = (eye(n_states) - K_state_ * H_state) * P_state_pred_damage;

    % Parameter prediction
    params_pred = params_est(:, k-1);
    P_param_pred = P_param / lambda + Q_param;
    
    % Parameter update using acceleration
    Z_accel = a_meas(:, k);
    x_vec = state_est(1:2:10, k);
    v_vec = state_est(2:2:10, k);
    a_pred = M \ (F_ext(:, k) - C_est*v_vec - K_est*x_vec);
    
    % Sensitivity matrix H_param (∂a/∂θ)
    H_param = zeros(5, n_params);
    
    % Stiffness sensitivities
    for i = 1:5
        dK_dki = zeros(5,5);
        if i == 1
            dK_dki(1,1) = 1; %  dK_dki(1,2) = 0; dK_dki(2,1) = 0; dK_dki(2,2) = 0;
        % elseif i == 5
            % dK_dki(4,4) = 1; dK_dki(4,5) = -1; dK_dki(5,4) = -1; dK_dki(5,5) = 1;
        else
            dK_dki(i,i) = 1; dK_dki(i,i-1) = -1; dK_dki(i-1,i) = -1 ;dK_dki(i-1,i-1) = 1;
            % dK_dki(i,i-1) = -1; dK_dki(i-1,i) = -1;
            %  dK_dki(i+1,i+1) = 1;
        end
        H_param(:, i) = -M \ (dK_dki * x_vec);
    end
    
    %H_param(:,1:5) = [x_vec(1), x_vec(1)-x_vec(2)]


    % Damping sensitivities
    for i = 1:5
        dC_dci = zeros(5,5);
        if i == 1
            dC_dci(1,1) = 1; %dC_dci(1,2) = -1; dC_dci(2,1) = -1; dC_dci(2,2) = 1;
        % elseif i == 5
        %     dC_dci(4,4) = 1; dC_dci(4,5) = -1; dC_dci(5,4) = -1; dC_dci(5,5) = 1;
        else
            dC_dci(i,i) = 1; dC_dci(i,i-1) = -1; dC_dci(i-1,i) = -1;dC_dci(i-1,i-1) = 1;
            % dC_dci(i,i-1) = -1; dC_dci(i-1,i) = -1;
            %  dC_dci(i+1,i+1) = 1;
        end
        H_param(:, i+5) = -M \ (dC_dci * v_vec);
    end
    
    innovation_param = Z_accel - a_pred;
    R_accel = diag([noise_acceleration^2 * mean(a_true(1,:).^2), ...
                    noise_acceleration^2 * mean(a_true(2,:).^2), ...
                    noise_acceleration^2 * mean(a_true(3,:).^2), ...
                    noise_acceleration^2 * mean(a_true(4,:).^2), ...
                    noise_acceleration^2 * mean(a_true(5,:).^2)]);
    
    S_param = H_param * P_param_pred * H_param' + R_accel;
    S_param = S_param + 1e-6 * eye(5);
    K_param = P_param_pred * H_param' / S_param;
    
    params_update = params_pred + K_param * innovation_param;
    
    % Apply constraints
    params_update(1:5) = max(params_update(1:5), 100);
    params_update(1:5) = min(params_update(1:5), 2000);
    params_update(6:10) = max(params_update(6:10), 0.1);
    params_update(6:10) = min(params_update(6:10), 20);
    
    params_est(:, k) = params_update;
    P_param = (eye(n_params) - K_param * H_param) * P_param_pred;
    
    % Store stiffness estimates
    k_est(:, k) = params_est(1:5, k);
    
    % ====================================================================
    % APPROACH 2: Damage factor as parameter (simplified)
    % ====================================================================


      % Damage factors evolve slowly (random walk)
    damage_pred = damage_est(:, k-1);

    P_param_damage_pred = P_param_damage + Q_param_damage;

    % Getting x, v vect 
    x_vec_ = state_est(1:2:10, k);
    v_vec_ = state_est(2:2:10, k);
    a_pred_ = M \ (F_ext(:, k) - C_est*v_vec_ - k_est_damage*x_vec_);

    % measurement  model for param fd
    innovation_param_damage = Z_accel - a_pred_;
    
     H_param_damage = [k1*x_vec_(1), k2*(x_vec_(1)-x_vec_(2)),0,0,0;
                        0, -k2*(x_vec_(1)-x_vec_(2)), k3*(x_vec_(2)-x_vec_(3)), 0, 0;
                        0, 0,-k3*(x_vec_(2)-x_vec_(3)), k4*(x_vec_(3)-x_vec_(4)), 0;
                         0, 0, 0, -k4*(x_vec_(3)-x_vec_(4)), k5*(x_vec_(4)-x_vec_(5));
                         0, 0, 0, 0, -k5*(x_vec_(4)-x_vec_(5))];


    % Param fd updation
    
    S_param_ = H_param_damage * P_param_damage_pred * H_param_damage' + R_accel;
    S_param_ = S_param_ + 1e-16 * eye(5);
    K_param_ = P_param_damage_pred * H_param_damage' / S_param_;
    
    damage_update = damage_pred + K_param_ * innovation_param_damage;

  
    % innovation_param_damage = 


    % Constrain damage factors to [0, 0.8]
    damage_update = max(damage_update, 0);
    damage_update = min(damage_update, 0.8);

    damage_est(:, k) = damage_update;

    % Calculate stiffness from damage factors
    k_est_damage_(1,k) = k1 * (1 - damage_est(1,k));
    k_est_damage_(2,k) = k2 * (1 - damage_est(2,k));
    k_est_damage_(3,k) = k3 * (1 - damage_est(3,k));
    k_est_damage_(4,k) = k4 * (1 - damage_est(4,k));
    k_est_damage_(5,k) = k5 * (1 - damage_est(5,k));
end

fprintf('DEKF completed.\n');

%% ========================================================================
% HELPER FUNCTIONS
% ========================================================================
function K = build_stiffness_matrix(k)
    % Build stiffness matrix for 5-DOF system
    K = zeros(5,5);
    K(1,1) = k(1) + k(2);
    K(1,2) = -k(2);
    K(2,1) = -k(2);
    K(2,2) = k(2) + k(3);
    K(2,3) = -k(3);
    K(3,2) = -k(3);
    K(3,3) = k(3) + k(4);
    K(3,4) = -k(4);
    K(4,3) = -k(4);
    K(4,4) = k(4) + k(5);
    K(4,5) = -k(5);
    K(5,4) = -k(5);
    K(5,5) = k(5);
end

function C = build_damping_matrix(c)
    % Build damping matrix for 5-DOF system
    C = zeros(5,5);
    C(1,1) = c(1) + c(2);
    C(1,2) = -c(2);
    C(2,1) = -c(2);
    C(2,2) = c(2) + c(3);
    C(2,3) = -c(3);
    C(3,2) = -c(3);
    C(3,3) = c(3) + c(4);
    C(3,4) = -c(4);
    C(4,3) = -c(4);
    C(4,4) = c(4) + c(5);
    C(4,5) = -c(5);
    C(5,4) = -c(5);
    C(5,5) = c(5);
end

%% ========================================================================
% VISUALIZATION - STIFFNESS TRACKING
% ========================================================================
fprintf('\nGenerating plots...\n');

% True stiffness values
k_true_all = [k1_true; k2_true; k3_true; k4_true; k5_true];

% Figure 1: All spring stiffness tracking (ki as parameters)
figure('Position', [50, 50, 1400, 1000]);

spring_names = {'k₁', 'k₂', 'k₃', 'k₄', 'k₅'};
initial_vals = [k1, k2, k3, k4, k5];

for i = 1:5
    subplot(5, 2, (i-1)*2+1);
    plot(time, k_true_all(i,:), 'k-', 'LineWidth', 3, 'DisplayName', sprintf('True %s', spring_names{i}));
    hold on;
    plot(time, k_est(i,:), 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Estimated %s', spring_names{i}));
    plot(time, initial_vals(i)*ones(size(time)), 'b:', 'LineWidth', 1.5, 'DisplayName', 'Initial');
    xlabel('Time (s)', 'FontSize', 11);
    ylabel(sprintf('%s (N/m)', spring_names{i}), 'FontSize', 11);
    title(sprintf('Spring %d Stiffness Tracking', i), 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    xlim([0, T]);
    ylim([0.5*initial_vals(i), 1.1*initial_vals(i)]);
    
    subplot(5, 2, (i-1)*2+2);
    error_k = abs(k_est(i,:) - k_true_all(i,:)) ./ k_true_all(i,:) * 100;
    plot(time, error_k, 'm-', 'LineWidth', 2);
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Estimation Error (%)', 'FontSize', 11);
    title(sprintf('Spring %d Estimation Error', i), 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    xlim([0, T]);
    ylim([0, 30]);
end

sgtitle('5-DOF System: Spring Stiffness Tracking (ki as parameters)', 'FontSize', 16, 'FontWeight', 'bold');

% Figure 2: Damage factor comparison
figure('Position', [50, 50, 1200, 800]);

% True damage factors
damage_true_all = [damage_k1; damage_k2; damage_k3; damage_k4; damage_k5];

for i = 1:5
    subplot(5, 2, (i-1)*2+1);
    plot(time, damage_true_all(i,:)*100, 'k-', 'LineWidth', 2, 'DisplayName', 'True Damage');
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Damage (%)', 'FontSize', 11);
    title(sprintf('Spring %d: True Damage', i), 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    xlim([0, T]);
    ylim([0, 70]);
    
    subplot(5, 2, (i-1)*2+2);
    est_damage_ki = 1 - k_est(i,:) / initial_vals(i);
    est_damage_fd = damage_est(i,:);
    
    plot(time, est_damage_ki*100, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Estimated from ki');
    hold on;
    plot(time, est_damage_fd*100, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Estimated damage factor');
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Damage (%)', 'FontSize', 11);
    title(sprintf('Spring %d: Estimated Damage', i), 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    xlim([0, T]);
    ylim([0, 70]);
end

sgtitle('5-DOF System: Damage Factor Comparison', 'FontSize', 16, 'FontWeight', 'bold');

% Figure 3: All stiffness in one plot
figure('Position', [50, 50, 1200, 500]);

colors = {'r', 'g', 'b', 'm', 'c'};
for i = 1:5
    plot(time, k_true_all(i,:), 'k-', 'LineWidth', 2);
    hold on;
    plot(time, k_est(i,:), '--', 'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Estimated %s', spring_names{i}));
    plot(time, k_est_damage_(i,:), ':', 'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Estimated %s', spring_names{i}));
end

xlabel('Time (s)', 'FontSize', 12);
ylabel('Spring Stiffness (N/m)', 'FontSize', 12);
title('5-DOF System: Simultaneous Damage Tracking', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, T]);

% Add damage percentage annotations
for i = 1:5
    text(2, initial_vals(i)*0.85, sprintf('%s: %.0f%% damage', spring_names{i}, max(damage_true_all(i,:))*100), ...
         'FontSize', 10, 'Color', colors{i}, 'BackgroundColor', 'white');
end

% Figure 4: State estimation performance
figure('Position', [50, 50, 1200, 600]);

subplot(2,1,1);
for i = 1:5
    plot(time, x_true(i,:), '-', 'LineWidth', 1, 'Color', colors{i}, ...
        'DisplayName', sprintf('True x_%d', i));
    hold on;
    plot(time, state_est(i*2-1,:), '--', 'LineWidth', 0.5, 'Color', colors{i}, ...
        'DisplayName', sprintf('Estimated x_%d', i));
end
xlabel('Time (s)', 'FontSize', 12);
ylabel('Displacement (m)', 'FontSize', 12);
title('Position Estimation', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, T]);

subplot(2,1,2);
for i = 1:5
    error_x = abs(state_est(i*2-1,:) - x_true(i,:));
    plot(time, error_x, '-', 'Color', colors{i}, 'LineWidth', 1, ...
        'DisplayName', sprintf('Error x_%d', i));
    hold on;
end
xlabel('Time (s)', 'FontSize', 12);
ylabel('Estimation Error (m)', 'FontSize', 12);
title('Position Estimation Errors', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, T]);

% Figure 5: Performance metrics
figure('Position', [50, 50, 1000, 400]);

% Calculate final errors and RMSE
final_errors = zeros(5,1);
rmse_values = zeros(5,1);

for i = 1:5
    final_errors(i) = abs(k_est(i,end) - k_true_all(i,end)) / k_true_all(i,end) * 100;
    rmse_values(i) = sqrt(mean((k_est(i,:) - k_true_all(i,:)).^2));
end

subplot(1, 2, 1);
bar(1:5, final_errors);
set(gca, 'XTick', 1:5, 'XTickLabel', spring_names);
ylabel('Final Estimation Error (%)', 'FontSize', 12);
title('Final Parameter Estimation Error', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ylim([0, 20]);

for i = 1:5
    text(i, final_errors(i)+0.5, sprintf('%.1f%%', final_errors(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end

subplot(1, 2, 2);
bar(1:5, rmse_values);
set(gca, 'XTick', 1:5, 'XTickLabel', spring_names);
ylabel('RMSE (N/m)', 'FontSize', 12);
title('Root Mean Square Error', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

for i = 1:5
    text(i, rmse_values(i)+5, sprintf('%.1f', rmse_values(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end



% figure 6 ( True vs actual in same plot )

figure('Position', [50, 50, 1200, 500]);

colors = {'r', 'g', 'b', 'm', 'c'};
for i = 1:5
    plot(time, damage_true_all(i,:),'-', 'Color', colors{i},'LineWidth', 2);
    hold on;
    plot(time, damage_est(i,:), '--', 'Color', colors{i}, 'LineWidth', 1, ...
        'DisplayName', sprintf('Estimated %s', spring_names{i}));
end

xlabel('Time (s)', 'FontSize', 12);
ylabel('Spring Stiffness (N/m)', 'FontSize', 12);
title('5-DOF System: Simultaneous Damage Tracking', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, T]);

%% ========================================================================
% PERFORMANCE ANALYSIS
% ========================================================================
fprintf('\n================================================================\n');
fprintf('PERFORMANCE ANALYSIS - 5-DOF SYSTEM\n');
fprintf('================================================================\n');

% RMSE calculations
for i = 1:5
    k_rmse = sqrt(mean((k_est(i,:) - k_true_all(i,:)).^2));
    fprintf('Spring %d (%s): RMSE = %.1f N/m (%.1f%% of initial value)\n', ...
            i, spring_names{i}, k_rmse, k_rmse/initial_vals(i)*100);
end

fprintf('\nFinal Stiffness Values (at t = %.1f s):\n', T);
fprintf('Spring\tTrue\t\tEstimated\tError\t\tDamage Detected\n');
for i = 1:5
    k_true_end = k_true_all(i,end);
    k_est_end = k_est(i,end);
    error = abs(k_est_end - k_true_end) / k_true_end * 100;
    damage_true = max(damage_true_all(i,:)) * 100;
    damage_detected = (1 - k_est_end/initial_vals(i)) * 100;
    
    fprintf('%d\t%.1f\t\t%.1f\t\t%.1f%%\t\t%.0f%% (True: %.0f%%)\n', ...
            i, k_true_end, k_est_end, error, damage_detected, damage_true);
end

% State estimation performance
state_rmse = sqrt(mean((state_est - state_true).^2, 2));
fprintf('\nState Estimation RMSE:\n');
for i = 1:5
    fprintf('  x%d: %.4f m\tv%d: %.4f m/s\n', i, state_rmse(i*2-1), i, state_rmse(i*2));
end

% Detection timing analysis
threshold = 0.95;  % 95% of initial value
true_start_times = [1, 5, 2, 6, 4];  % Approximate start times

fprintf('\nDamage Detection Timing (when stiffness drops below 95%% of initial):\n');
fprintf('Spring\tTrue Start\tDetection\tDelay\n');
for i = 1:5
    idx = find(k_est(i,:) < threshold * initial_vals(i), 1);
    if ~isempty(idx)
        detection_time = time(idx);
        delay = detection_time - true_start_times(i);
        fprintf('%s\t%.1f s\t\t%.1f s\t\t%.1f s\n', spring_names{i}, true_start_times(i), detection_time, delay);
    else
        fprintf('%s\t%.1f s\t\tNot detected\t-\n', spring_names{i}, true_start_times(i));
    end
end

fprintf('\n================================================================\n');
fprintf('Simulation Complete\n');
fprintf('Total time: %.1f s, Time step: %.4f s\n', T, dt);
fprintf('5-DOF System with %d springs tracked\n', 5);
fprintf('================================================================\n');
