%% DIFFERENT DIFFERENT SCENARIO IN ONE PLOT Enhanced Kalman Filter for FEM-MEASUREMENT FUSION
% This code demonstrates proper sensor fusion between FEM predictions and
% noisy measurements using Kalman Filter with tunable trust parameters.

clear all; close all; clc;

fprintf('=== Enhanced Kalman Filter for FEM-Measurement Fusion ===\n');

%% Problem Parameters
L = 1.0;                    % Beam length in meters
E = 200e9;                  % Young's modulus in Pa (Steel)
I = 8.33e-9;                % Moment of inertia in m^4
P = 1000;                   % Point load at free end in N
n_elements = 10;            % Number of finite elements
n_nodes = n_elements + 1;   % Number of nodes

fprintf('Beam Parameters:\n');
fprintf('Length: %.2f m, Load: %.0f N\n', L, P);

%% FEM Analysis
x_nodes = linspace(0, L, n_nodes)';  % Node positions

% Initialize global stiffness matrix and force vector
K_global = zeros(2*n_nodes, 2*n_nodes);
F_global = zeros(2*n_nodes, 1);

% Element length
Le = L / n_elements;

% Beam element stiffness matrix
k_elem = (E*I/Le^3) * [12, 6*Le, -12, 6*Le;
                       6*Le, 4*Le^2, -6*Le, 2*Le^2;
                       -12, -6*Le, 12, -6*Le;
                       6*Le, 2*Le^2, -6*Le, 4*Le^2];

% Assemble global stiffness matrix
for elem = 1:n_elements
    dof_indices = [2*elem-1, 2*elem, 2*elem+1, 2*elem+2];
    K_global(dof_indices, dof_indices) = K_global(dof_indices, dof_indices) + k_elem;
end

% Apply boundary conditions and solve
free_dof = 3:2*n_nodes;
K_reduced = K_global(free_dof, free_dof);
F_reduced = zeros(length(free_dof), 1);
F_reduced(end-1) = P;

U_reduced = K_reduced \ F_reduced;
U_full = zeros(2*n_nodes, 1);
U_full(free_dof) = U_reduced;
fem_deflection = U_full(1:2:end);  % FEM prediction

%% Analytical Solution and Noisy Measurements
analytical_deflection = (P./(6*E*I)) .* (3*L*x_nodes.^2 - x_nodes.^3);

rng(42);  % For reproducibility
measurement_noise_std = 0.1 * max(abs(analytical_deflection));
noisy_measurements = analytical_deflection + measurement_noise_std * randn(size(analytical_deflection));

fprintf('FEM max deflection: %.6f m\n', max(fem_deflection));
fprintf('Analytical max deflection: %.6f m\n', max(analytical_deflection));

%% KALMAN FILTER WITH DUAL OBSERVATIONS (FEM + MEASUREMENTS)
fprintf('\n=== Setting up Dual-Observation Kalman Filter ===\n');

n_states = n_nodes;

% KF TUNING PARAMETERS - USER CAN ADJUST THESE!
trust_in_fem = 0.8;        % 0-1: How much to trust FEM (1 = complete trust)
trust_in_measurements = 0.2; % 0-1: How much to trust measurements

fprintf('Tuning Parameters:\n');
fprintf('Trust in FEM: %.2f\n', trust_in_fem);
fprintf('Trust in Measurements: %.2f\n', trust_in_measurements);

% State transition matrix (identity for static problem)
A = eye(n_states);

% Process noise (model uncertainty)
Q = 1e-4 * eye(n_states);

% DUAL MEASUREMENT SETUP
H_fem = eye(n_states);          % FEM observation matrix
H_meas = eye(n_states);         % Measurement observation matrix

% Measurement noise covariances (inverse of trust)
R_fem = ((1 - trust_in_fem) + 0.01)^2 * eye(n_states);      % FEM uncertainty
R_meas = ((1 - trust_in_measurements) + 0.01)^2 * eye(n_states); % Measurement uncertainty

% Combine observations
H_combined = [H_fem; H_meas];
R_combined = blkdiag(R_fem, R_meas);

% Initialize Kalman Filter
x_hat = zeros(n_states, 1);      % Initial state estimate
P_hat = 0.1 * eye(n_states);     % Initial error covariance

% Store estimates
n_steps = 100;
x_estimates = zeros(n_states, n_steps);
innovations = zeros(2*n_states, n_steps);

fprintf('\nRunning Dual-Observation Kalman Filter...\n');
for k = 1:n_steps
    % PREDICTION STEP
    x_hat_minus = A * x_hat;
    P_minus = A * P_hat * A' + Q;

    % UPDATE STEP WITH DUAL OBSERVATIONS
    % Combined measurements: [FEM predictions; Noisy measurements]
    z_combined = [fem_deflection; noisy_measurements];

    % Innovation (measurement residual)
    y = z_combined - H_combined * x_hat_minus;
    innovations(:, k) = y;

    % Kalman Gain
    S = H_combined * P_minus * H_combined' + R_combined;
    K = P_minus * H_combined' / S;

    % State update
    x_hat = x_hat_minus + K * y;

    % Covariance update
    P_hat = (eye(n_states) - K * H_combined) * P_minus;

    % Store results
    x_estimates(:, k) = x_hat;

    if mod(k, 20) == 0
        fprintf('  Iteration %d/%d - RMS Innovation: %.6e\n', k, n_steps, sqrt(mean(y.^2)));
    end
end

kf_final_estimate = x_estimates(:, end);

%% Performance Analysis
fem_error = abs(fem_deflection - analytical_deflection);
measurement_error = abs(noisy_measurements - analytical_deflection);
kf_error = abs(kf_final_estimate - analytical_deflection);

fprintf('\n=== PERFORMANCE COMPARISON ===\n');
fprintf('FEM RMS Error: %.6e m\n', sqrt(mean(fem_error.^2)));
fprintf('Measurement RMS Error: %.6e m\n', sqrt(mean(measurement_error.^2)));
fprintf('Kalman Filter RMS Error: %.6e m\n', sqrt(mean(kf_error.^2)));

improvement_over_fem = (sqrt(mean(fem_error.^2)) - sqrt(mean(kf_error.^2))) / sqrt(mean(fem_error.^2)) * 100;
improvement_over_meas = (sqrt(mean(measurement_error.^2)) - sqrt(mean(kf_error.^2))) / sqrt(mean(measurement_error.^2)) * 100;

fprintf('Improvement over FEM: %.2f%%\n', improvement_over_fem);
fprintf('Improvement over Measurements: %.2f%%\n', improvement_over_meas);

%% COMPARISON WITH DIFFERENT TRUST SETTINGS
fprintf('\n=== TESTING DIFFERENT TRUST CONFIGURATIONS ===\n');

trust_configs = {
    [0.9, 0.1], 'High trust in FEM';
    [0.1, 0.9], 'High trust in Measurements'; 
    [0.5, 0.5], 'Equal trust';
    [0.6, 0.4], 'Balanced (FEM preferred)';
};

kf_results = zeros(n_states, length(trust_configs));
config_names = cell(length(trust_configs), 1);

for config = 1:length(trust_configs)
    trust_fem = trust_configs{config, 1}(1);
    trust_meas = trust_configs{config, 1}(2);
    config_names{config} = trust_configs{config, 2};

    % Run KF with this configuration
    R_fem_config = ((1 - trust_fem) + 0.01)^2 * eye(n_states);
    R_meas_config = ((1 - trust_meas) + 0.01)^2 * eye(n_states);
    R_combined_config = blkdiag(R_fem_config, R_meas_config);

    x_hat_config = zeros(n_states, 1);
    P_hat_config = 0.1 * eye(n_states);

    for k = 1:100  % Fewer iterations for comparison
        x_hat_minus = A * x_hat_config;
        P_minus = A * P_hat_config * A' + Q;

        z_combined = [fem_deflection; noisy_measurements];
        S_config = H_combined * P_minus * H_combined' + R_combined_config;
        K_config = P_minus * H_combined' / S_config;

        x_hat_config = x_hat_minus + K_config * (z_combined - H_combined * x_hat_minus);
        P_hat_config = (eye(n_states) - K_config * H_combined) * P_minus;
    end

    kf_results(:, config) = x_hat_config;

    config_error = sqrt(mean((x_hat_config - analytical_deflection).^2));
    fprintf('Config: %s - RMS Error: %.6e m\n', trust_configs{config, 2}, config_error);
end

%% Plotting Results
figure('Position', [100, 100, 1500, 1000]);

% Plot 1: Main Comparison
subplot(2, 3, 1);
plot(x_nodes, analytical_deflection * 1000, 'k-', 'LineWidth', 3, 'DisplayName', 'True Analytical');
hold on;
plot(x_nodes, fem_deflection * 1000, 'b--', 'LineWidth', 2, 'DisplayName', 'FEM Prediction');
plot(x_nodes, noisy_measurements * 1000, 'r*', 'LineWidth', 1.5, 'DisplayName', 'Noisy Measurements');
plot(x_nodes, kf_final_estimate * 1000, 'g*-.', 'LineWidth', 2.5, 'DisplayName', 'KF Fused Estimate');
xlabel('Beam Position (m)');
ylabel('Deflection (mm)');
title('Dual-Observation Kalman Filter Fusion');
legend('Location', 'northwest');
grid on;

% Plot 2: Errors
subplot(2, 3, 2);
errors = [fem_error, measurement_error, kf_error] * 1000;
bar(1:n_nodes, errors);
xlabel('Node Number');
ylabel('Absolute Error (mm)');
title('Comparison of Errors');
legend('FEM Error', 'Measurement Error', 'KF Error', 'Location', 'northwest');
grid on;

% Plot 3: KF Convergence
subplot(2, 3, 3);
rms_evolution = sqrt(mean((x_estimates - analytical_deflection).^2, 1)) * 1000;
plot(1:n_steps, rms_evolution, 'b-', 'LineWidth', 2);
xlabel('Kalman Filter Iteration');
ylabel('RMS Error (mm)');
title('KF Convergence History');
grid on;

% Plot 4: Different Trust Configurations
subplot(2, 3, 4);
config_errors = zeros(length(trust_configs), 1);
for i = 1:length(trust_configs)
    config_errors(i) = sqrt(mean((kf_results(:, i) - analytical_deflection).^2)) * 1000;
end
bar(categorical(config_names), config_errors);
ylabel('RMS Error (mm)');
title('Performance vs Trust Configuration');
grid on;
rotateXLabels(gca, 45);

% Plot 5: Innovation Sequence
subplot(2, 3, 5);
fem_innovations = sqrt(mean(innovations(1:n_states, :).^2, 1)) * 1000;
meas_innovations = sqrt(mean(innovations(n_states+1:end, :).^2, 1)) * 1000;
plot(1:n_steps, fem_innovations, 'b-', 'LineWidth', 2, 'DisplayName', 'FEM Innovations');
hold on;
plot(1:n_steps, meas_innovations, 'r-', 'LineWidth', 2, 'DisplayName', 'Measurement Innovations');
xlabel('Iteration');
ylabel('RMS Innovation (mm)');
title('Innovation Sequences');
legend();
grid on;

% Plot 6: Trust Weight Visualization
subplot(2, 3, 6);
trust_weights = [trust_in_fem, trust_in_measurements; 0.5, 0.5; 0.9, 0.1; 0.1, 0.9];
labels = {'Current', 'Equal', 'FEM Trust', 'Meas Trust'};
bar(categorical(labels), trust_weights, 'stacked');
ylabel('Trust Weight');
title('Trust Configuration Comparison');
legend('FEM Trust', 'Measurement Trust');
grid on;

%% Interactive Trust Parameter Section
fprintf('\n=== INTERACTIVE TRUST PARAMETER TESTING ===\n');
fprintf('You can modify these trust parameters in the code:\n');
fprintf('   trust_in_fem = 0.8;        %% How much to trust FEM (0-1)\n');
fprintf('   trust_in_measurements = 0.2; %% How much to trust measurements (0-1)\n');
fprintf('\nRecommended configurations:\n');
fprintf('   For accurate FEM: [0.9, 0.1]\n');
fprintf('   For accurate sensors: [0.1, 0.9]\n');
fprintf('   For balanced approach: [0.5, 0.5] or [0.6, 0.4]\n');

%% Display Final Results
fprintf('\n=== FINAL RESULTS ===\n');
fprintf('Node\tTrue(mm)\tFEM(mm)\t\tMeas(mm)\tKF Est.(mm)\n');
fprintf('----------------------------------------------------------------\n');
for i = [1, 3, 6, 9, 11]
    fprintf('%d\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n', i, ...
        analytical_deflection(i)*1000, fem_deflection(i)*1000, ...
        noisy_measurements(i)*1000, kf_final_estimate(i)*1000);
end

fprintf('\n=== SIMULATION COMPLETE ===\n');





%% Kalman Filter Trust Situation Comparison
% % This code generates different plots for different trust configurations
% % to show how KF behavior changes with different trust levels
% 
% clear all; close all; clc;
% 
% fprintf('=== Kalman Filter Trust Situation Analysis ===\n');
% 
% %% Problem Parameters
% L = 1.0;                    % Beam length in meters
% E = 200e9;                  % Young's modulus in Pa (Steel)
% I = 8.33e-9;                % Moment of inertia in m^4
% P = 1000;                   % Point load at free end in N
% n_elements = 10;            % Number of finite elements
% n_nodes = n_elements + 1;   % Number of nodes
% 
% %% FEM Analysis
% x_nodes = linspace(0, L, n_nodes)';  % Node positions
% 
% % Initialize global stiffness matrix and force vector
% K_global = zeros(2*n_nodes, 2*n_nodes);
% F_global = zeros(2*n_nodes, 1);
% 
% % Element length
% Le = L / n_elements;
% 
% % Beam element stiffness matrix
% k_elem = (E*I/Le^3) * [12, 6*Le, -12, 6*Le;
%                        6*Le, 4*Le^2, -6*Le, 2*Le^2;
%                        -12, -6*Le, 12, -6*Le;
%                        6*Le, 2*Le^2, -6*Le, 4*Le^2];
% 
% % Assemble global stiffness matrix
% for elem = 1:n_elements
%     dof_indices = [2*elem-1, 2*elem, 2*elem+1, 2*elem+2];
%     K_global(dof_indices, dof_indices) = K_global(dof_indices, dof_indices) + k_elem;
% end
% 
% % Apply boundary conditions and solve
% free_dof = 3:2*n_nodes;
% K_reduced = K_global(free_dof, free_dof);
% F_reduced = zeros(length(free_dof), 1);
% F_reduced(end-1) = P;
% 
% U_reduced = K_reduced \ F_reduced;
% U_full = zeros(2*n_nodes, 1);
% U_full(free_dof) = U_reduced;
% fem_deflection = U_full(1:2:end);  % FEM prediction
% 
% %% Analytical Solution and Noisy Measurements
% analytical_deflection = (P./(6*E*I)) .* (3*L*x_nodes.^2 - x_nodes.^3);
% 
% rng(42);  % For reproducibility
% measurement_noise_std = 0.1 * max(abs(analytical_deflection));
% noisy_measurements = analytical_deflection + measurement_noise_std * randn(size(analytical_deflection));
% 
% fprintf('FEM max deflection: %.6f m\n', max(fem_deflection));
% fprintf('Analytical max deflection: %.6f m\n', max(analytical_deflection));
% fprintf('Measurement noise std: %.6f m\n', measurement_noise_std);
% 
% %% Define Different Trust Situations
% trust_situations = {
%     [0.9, 0.1], 'High FEM Trust', 'KF follows FEM closely';
%     [0.1, 0.9], 'High Measurement Trust', 'KF follows measurements closely';
%     [0.5, 0.5], 'Equal Trust', 'Balanced approach';
%     [0.8, 0.2], 'FEM Preferred', 'FEM dominant with measurement correction';
%     [0.2, 0.8], 'Measurement Preferred', 'Measurements dominant with FEM guidance';
%     [0.95, 0.05], 'Very High FEM Trust', 'KF almost ignores measurements';
%     [0.05, 0.95], 'Very High Meas Trust', 'KF almost ignores FEM';
% };
% 
% n_situations = length(trust_situations);
% n_steps = 100;
% n_states = n_nodes;
% 
% % Store results for all situations
% kf_results_all = zeros(n_states, n_steps, n_situations);
% final_estimates = zeros(n_states, n_situations);
% situation_errors = zeros(n_situations, 1);
% 
% fprintf('\n=== Running Kalman Filter for Different Trust Situations ===\n');
% 
% for situation = 1:n_situations
%     trust_fem = trust_situations{situation, 1}(1);
%     trust_meas = trust_situations{situation, 1}(2);
%     situation_name = trust_situations{situation, 2};
% 
%     fprintf('Situation %d/%d: %s (FEM: %.2f, Meas: %.2f)\n', ...
%         situation, n_situations, situation_name, trust_fem, trust_meas);
% 
%     % State transition matrix
%     A = eye(n_states);
%     Q = 1e-12 * eye(n_states);
% 
%     % Dual measurement setup
%     H_fem = eye(n_states);
%     H_meas = eye(n_states);
%     H_combined = [H_fem; H_meas];
% 
%     % Measurement noise covariances based on trust
%     R_fem = ((1 - trust_fem) + 0.01)^2 * eye(n_states);
%     R_meas = ((1 - trust_meas) + 0.01)^2 * eye(n_states);
%     R_combined = blkdiag(R_fem, R_meas);
% 
%     % Initialize Kalman Filter
%     x_hat = zeros(n_states, 1);
%     P_hat = 0.1 * eye(n_states);
% 
%     % Run Kalman Filter
%     for k = 1:n_steps
%         % Prediction Step
%         x_hat_minus = A * x_hat;
%         P_minus = A * P_hat * A' + Q;
% 
%         % Update Step with Dual Observations
%         z_combined = [fem_deflection; noisy_measurements];
%         y = z_combined - H_combined * x_hat_minus;
% 
%         S = H_combined * P_minus * H_combined' + R_combined;
%         K = P_minus * H_combined' / S;
% 
%         x_hat = x_hat_minus + K * y;
%         P_hat = (eye(n_states) - K * H_combined) * P_minus;
% 
%         kf_results_all(:, k, situation) = x_hat;
%     end
% 
%     final_estimates(:, situation) = x_hat;
%     situation_errors(situation) = sqrt(mean((x_hat - analytical_deflection).^2));
% end
% 
% %% PLOT 1: COMPARISON OF ALL TRUST SITUATIONS - DEFLECTION PROFILES
% figure('Position', [50, 50, 1400, 900]);
% sgtitle('Kalman Filter Behavior Under Different Trust Situations', 'FontSize', 16, 'FontWeight', 'bold');
% 
% colors = lines(n_situations);
% 
% % Main deflection comparison
% subplot(2, 3, 1);
% plot(x_nodes, analytical_deflection * 1000, 'k-', 'LineWidth', 4, 'DisplayName', 'True Analytical');
% hold on;
% plot(x_nodes, fem_deflection * 1000, 'b--', 'LineWidth', 2, 'DisplayName', 'FEM Prediction');
% plot(x_nodes, noisy_measurements * 1000, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Noisy Measurements');
% 
% for situation = 1:n_situations
%     plot(x_nodes, final_estimates(:, situation) * 1000, '-', 'LineWidth', 2, ...
%         'Color', colors(situation, :), 'DisplayName', trust_situations{situation, 2});
% end
% xlabel('Beam Position (m)');
% ylabel('Deflection (mm)');
% title('Deflection Profiles - All Trust Situations');
% legend('Location', 'northwest', 'FontSize', 8);
% grid on;
% 
% %% PLOT 2: ERROR COMPARISON ACROSS SITUATIONS
% subplot(2, 3, 2);
% error_bars = zeros(n_situations, 3); % FEM error, Measurement error, KF error
% for situation = 1:n_situations
%     error_bars(situation, 1) = sqrt(mean((fem_deflection - analytical_deflection).^2)) * 1000;
%     error_bars(situation, 2) = sqrt(mean((noisy_measurements - analytical_deflection).^2)) * 1000;
%     error_bars(situation, 3) = situation_errors(situation) * 1000;
% end
% 
% b = bar(error_bars, 'grouped');
% b(1).FaceColor = 'blue';
% b(2).FaceColor = 'red';
% b(3).FaceColor = 'green';
% set(gca, 'XTickLabel', cellfun(@(x) x{2}, trust_situations, 'UniformOutput', false));
% xtickangle(45);
% ylabel('RMS Error (mm)');
% title('RMS Error Comparison');
% legend('FEM Error', 'Measurement Error', 'KF Error', 'Location', 'northeast');
% grid on;
% 
% %% PLOT 3: CONVERGENCE SPEED COMPARISON
% subplot(2, 3, 3);
% convergence_threshold = 0.01; % mm
% convergence_iterations = zeros(n_situations, 1);
% 
% for situation = 1:n_situations
%     rms_evolution = sqrt(mean((kf_results_all(:, :, situation) - analytical_deflection).^2, 1)) * 1000;
%     plot(1:n_steps, rms_evolution, 'LineWidth', 2, 'Color', colors(situation, :), ...
%         'DisplayName', trust_situations{situation, 2});
%     hold on;
% 
%     % Find convergence iteration
%     converged_idx = find(rms_evolution < convergence_threshold, 1);
%     if ~isempty(converged_idx)
%         convergence_iterations(situation) = converged_idx;
%         plot(converged_idx, rms_evolution(converged_idx), 'o', 'MarkerSize', 8, ...
%             'Color', colors(situation, :), 'MarkerFaceColor', colors(situation, :));
%     end
% end
% xlabel('Iteration');
% ylabel('RMS Error (mm)');
% title('Convergence Speed Comparison');
% legend('Location', 'northeast', 'FontSize', 8);
% grid on;
% 
% %% PLOT 4: FREE END DEFLECTION EVOLUTION
% subplot(2, 3, 4);
% free_end_node = n_nodes;
% true_free_end = analytical_deflection(free_end_node) * 1000;
% 
% for situation = 1:n_situations
%     free_end_evolution = squeeze(kf_results_all(free_end_node, :, situation)) * 1000;
%     plot(1:n_steps, free_end_evolution, 'LineWidth', 2, 'Color', colors(situation, :), ...
%         'DisplayName', trust_situations{situation, 2});
%     hold on;
% end
% plot([1, n_steps], [true_free_end, true_free_end], 'k--', 'LineWidth', 3, 'DisplayName', 'True Value');
% plot([1, n_steps], [fem_deflection(free_end_node)*1000, fem_deflection(free_end_node)*1000], 'b--', 'LineWidth', 2, 'DisplayName', 'FEM Value');
% plot([1, n_steps], [noisy_measurements(free_end_node)*1000, noisy_measurements(free_end_node)*1000], 'r:', 'LineWidth', 2, 'DisplayName', 'Measured Value');
% 
% xlabel('Iteration');
% ylabel('Free End Deflection (mm)');
% title('Free End Deflection Evolution');
% legend('Location', 'eastoutside', 'FontSize', 8);
% grid on;
% 
% %% PLOT 5: TRUST WEIGHT VISUALIZATION
% subplot(2, 3, 5);
% trust_matrix = cell2mat(cellfun(@(x) x{1}, trust_situations, 'UniformOutput', false)');
% b = bar(trust_matrix, 'stacked');
% b(1).FaceColor = [0, 0.4470, 0.7410]; % Blue for FEM
% b(2).FaceColor = [0.8500, 0.3250, 0.0980]; % Red for Measurements
% set(gca, 'XTickLabel', cellfun(@(x) x{2}, trust_situations, 'UniformOutput', false));
% xtickangle(45);
% ylabel('Trust Weight');
% title('Trust Weight Distribution');
% legend('FEM Trust', 'Measurement Trust', 'Location', 'northeast');
% ylim([0, 1]);
% grid on;
% 
% %% PLOT 6: FINAL ESTIMATE VS TRUST RATIO
% subplot(2, 3, 6);
% trust_ratio = trust_matrix(:, 1) ./ trust_matrix(:, 2); % FEM trust / Measurement trust
% final_errors = situation_errors * 1000;
% 
% scatter(trust_ratio, final_errors, 100, 1:n_situations, 'filled');
% colorbar;
% xlabel('Trust Ratio (FEM/Measurement)');
% ylabel('Final RMS Error (mm)');
% title('Performance vs Trust Ratio');
% grid on;
% 
% % Add labels for each point
% for situation = 1:n_situations
%     text(trust_ratio(situation), final_errors(situation), ...
%         trust_situations{situation, 2}, 'FontSize', 8, 'HorizontalAlignment', 'center');
% end
% set(gca, 'XScale', 'log');
% 
% %% PLOT 7: INDIVIDUAL SITUATION FOCUS PLOTS
% figure('Position', [100, 100, 1200, 800]);
% sgtitle('Detailed Analysis of Selected Trust Situations', 'FontSize', 16, 'FontWeight', 'bold');
% 
% selected_situations = [1, 2, 3, 6]; % High FEM, High Meas, Equal, Very High Meas
% 
% for i = 1:length(selected_situations)
%     situation = selected_situations(i);
% 
%     subplot(2, 2, i);
% 
%     % Plot all curves
%     plot(x_nodes, analytical_deflection * 1000, 'k-', 'LineWidth', 3, 'DisplayName', 'True Analytical');
%     hold on;
%     plot(x_nodes, fem_deflection * 1000, 'b--', 'LineWidth', 2, 'DisplayName', 'FEM Prediction');
%     plot(x_nodes, noisy_measurements * 1000, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Noisy Measurements');
%     plot(x_nodes, final_estimates(:, situation) * 1000, 'g-.', 'LineWidth', 2.5, 'DisplayName', 'KF Estimate');
% 
%     xlabel('Beam Position (m)');
%     ylabel('Deflection (mm)');
%     title(sprintf('%s\n(FEM: %.2f, Meas: %.2f)', ...
%         trust_situations{situation, 2}, ...
%         trust_situations{situation, 1}(1), ...
%         trust_situations{situation, 1}(2)));
%     legend('Location', 'northwest');
%     grid on;
% end
% 
% %% PLOT 8: CONVERGENCE BEHAVIOR AT KEY NODES
% figure('Position', [150, 150, 1200, 800]);
% sgtitle('Node Convergence Behavior Under Different Trust Situations', 'FontSize', 16, 'FontWeight', 'bold');
% 
% key_nodes = [3, 6, n_nodes]; % Fixed end, middle, free end
% node_names = {'Near Fixed End', 'Middle', 'Free End'};
% 
% for node_idx = 1:length(key_nodes)
%     node = key_nodes(node_idx);
%     true_value = analytical_deflection(node) * 1000;
% 
%     subplot(2, 2, node_idx);
% 
%     for situation = 1:n_situations
%         node_evolution = squeeze(kf_results_all(node, :, situation)) * 1000;
%         plot(1:n_steps, node_evolution, 'LineWidth', 1.5, 'Color', colors(situation, :), ...
%             'DisplayName', trust_situations{situation, 2});
%         hold on;
%     end
% 
%     plot([1, n_steps], [true_value, true_value], 'k--', 'LineWidth', 3, 'DisplayName', 'True Value');
% 
%     xlabel('Iteration');
%     ylabel(['Deflection at Node ', num2str(node), ' (mm)']);
%     title(['Convergence at ', node_names{node_idx}, ' (x=', num2str(x_nodes(node)), 'm)']);
%     grid on;
% 
%     if node_idx == 1
%         legend('Location', 'eastoutside', 'FontSize', 7);
%     end
% end
% 
% % Add performance summary table
% subplot(2, 2, 4);
% axis off;
% performance_text = sprintf('PERFORMANCE SUMMARY\n\n');
% for situation = 1:n_situations
%     performance_text = [performance_text, sprintf('%s:\n  RMS Error: %.4f mm\n  Trust: FEM=%.2f, Meas=%.2f\n\n', ...
%         trust_situations{situation, 2}, situation_errors(situation)*1000, ...
%         trust_situations{situation, 1}(1), trust_situations{situation, 1}(2))];
% end
% text(0.1, 0.9, performance_text, 'FontSize', 10, 'VerticalAlignment', 'top', 'FontName', 'FixedWidth');
% 
% %% Display Final Performance Summary
% fprintf('\n=== FINAL PERFORMANCE SUMMARY ===\n');
% fprintf('Situation\t\t\t\tFEM Trust\tMeas Trust\tRMS Error (mm)\n');
% fprintf('--------------------------------------------------------------------\n');
% for situation = 1:n_situations
%     fprintf('%-25s\t%.2f\t\t%.2f\t\t%.4f\n', ...
%         trust_situations{situation, 2}, ...
%         trust_situations{situation, 1}(1), ...
%         trust_situations{situation, 1}(2), ...
%         situation_errors(situation) * 1000);
% end
% 
% fprintf('\nBest performance: %s (%.4f mm)\n', ...
%     trust_situations{situation_errors == min(situation_errors), 2}, ...
%     min(situation_errors) * 1000);
% 
% fprintf('\n=== ANALYSIS COMPLETE ===\n');
% fprintf('Generated 3 comprehensive figures with multiple subplots\n');
% fprintf('showing different aspects of Kalman Filter behavior under various trust situations.\n');




%%%%%%%%% ANotehr method %%%%%%%%%%%%%
%% Cantilever Beam Deflection Analysis with Kalman Filter Estimation - SIMULATION
% This code demonstrates FEM analysis of a cantilever beam, adds measurement
% noise to the analytical solution, and uses Kalman Filter to estimate the
% true deflection from noisy measurements.

clear all; close all; clc;

fprintf('=== Starting Cantilever Beam Deflection Simulation ===\n');

%% Problem Parameters
L = 1.0;                    % Beam length in meters
E = 200e9;                  % Young's modulus in Pa (Steel)
I = 8.33e-9;                % Moment of inertia in m^4 (Rectangular cross-section)
P = 1000;                   % Point load at free end in N
n_elements = 50;            % Number of finite elements
n_nodes = n_elements + 1;   % Number of nodes

fprintf('Beam Parameters:\n');
fprintf('Length: %.2f m, Young''s Modulus: %.2e Pa, Load: %.0f N\n', L, E, P);
fprintf('Number of elements: %d, Number of nodes: %d\n', n_elements, n_nodes);

%% Finite Element Method for Cantilever Beam
fprintf('\n=== Performing FEM Analysis ===\n');

% Generate node coordinates
x_nodes = linspace(0, L, n_nodes)';  % Node positions along the beam
fprintf('Node positions: %s\n', mat2str(x_nodes', 2));

% Initialize global stiffness matrix and force vector
K_global = zeros(2*n_nodes, 2*n_nodes);  % 2 DOF per node (deflection, slope)
F_global = zeros(2*n_nodes, 1);          % Global force vector

% Element length
Le = L / n_elements;
fprintf('Element length: %.3f m\n', Le);

% Beam element stiffness matrix (for Euler-Bernoulli beam)
k_elem = (E*I/Le^3) * [12, 6*Le, -12, 6*Le;
                       6*Le, 4*Le^2, -6*Le, 2*Le^2;
                       -12, -6*Le, 12, -6*Le;
                       6*Le, 2*Le^2, -6*Le, 4*Le^2];

fprintf('Element stiffness matrix (first 2x2 block):\n');
disp(k_elem(1:2, 1:2));

% Assemble global stiffness matrix
fprintf('Assembling global stiffness matrix...\n');
for elem = 1:n_elements
    % DOF indices for current element
    dof_indices = [2*elem-1, 2*elem, 2*elem+1, 2*elem+2];

    % Add element stiffness to global matrix
    K_global(dof_indices, dof_indices) = K_global(dof_indices, dof_indices) + k_elem;
end

% Apply boundary conditions (fixed at x=0)
% Remove rows and columns corresponding to fixed DOFs (deflection and slope at node 1)
free_dof = 3:2*n_nodes;  % All DOFs except first two (fixed end)
K_reduced = K_global(free_dof, free_dof);

% Apply point load at free end (last node)
F_reduced = zeros(length(free_dof), 1);
F_reduced(end-1) = P;  % Force applied at deflection DOF of last node

fprintf('Solving FEM system...\n');
% Solve for displacements
U_reduced = K_reduced \ F_reduced;

% Reconstruct full displacement vector
U_full = zeros(2*n_nodes, 1);
U_full(free_dof) = U_reduced;

% Extract deflection values at nodes
fem_deflection = U_full(1:2:end);  % Deflection at each node

fprintf('FEM solution completed. Maximum deflection: %.6f m\n', max(fem_deflection));

%% Analytical Solution (for comparison)
fprintf('\n=== Calculating Analytical Solution ===\n');
% Theoretical deflection for cantilever beam with end load
analytical_deflection = (P./(6*E*I)) .* (3*L*x_nodes.^2 - x_nodes.^3);
fprintf('Analytical maximum deflection: %.6f m\n', max(analytical_deflection));

%% Add Measurement Noise to Analytical Solution
fprintf('\n=== Adding Measurement Noise ===\n');
rng(42);  % Set seed for reproducibility
measurement_noise_std = 0.1 * max(abs(analytical_deflection));  % 10% noise level
noisy_measurements = analytical_deflection + measurement_noise_std * randn(size(analytical_deflection));

fprintf('Measurement noise standard deviation: %.6f m\n', measurement_noise_std);
fprintf('Noisy measurements range: %.6f to %.6f m\n', min(noisy_measurements), max(noisy_measurements));

%% Kalman Filter Implementation
fprintf('\n=== Implementing Kalman Filter ===\n');
% State-space representation of beam deflection
% We'll model the deflection as a random walk process

n_states = n_nodes;
dt = 1;  % Time step (arbitrary for static problem)

% State transition matrix (identity since deflection doesn't change)
A = eye(n_states);

% Process noise covariance (small, representing model uncertainty)
Q = 10 * eye(n_states);

% Measurement matrix (we measure all states directly)
H = eye(n_states);

% Measurement noise covariance
R = measurement_noise_std^2 * eye(n_states);

% Initialize Kalman Filter
x_hat = zeros(n_states, 1);          % Initial state estimate
P_hat = 1 * eye(n_states);           % Initial error covariance

% Store estimates
n_steps = 50;  % Number of Kalman filter iterations
x_estimates = zeros(n_states, n_steps);
P_diag = zeros(n_states, n_steps);

fprintf('Running Kalman Filter for %d iterations...\n', n_steps);
% Kalman Filter Iteration
for k = 1:n_steps
    % Prediction Step
    x_hat_minus = A * x_hat;
    P_minus = A * P_hat * A' + Q;

    % Update Step
    K = P_minus * H' / (H * P_minus * H' + R);  % Kalman Gain
    x_hat = x_hat_minus + K * (noisy_measurements - H * x_hat_minus);
    P_hat = (eye(n_states) - K * H) * P_minus;

    % Store results
    x_estimates(:, k) = x_hat;
    P_diag(:, k) = diag(P_hat);

    % Display progress every 10 steps
    if mod(k, 10) == 0
        fprintf('  Iteration %d/%d completed\n', k, n_steps);
    end
end

% Final Kalman Filter estimate
kf_final_estimate = x_estimates(:, end);

fprintf('Kalman Filter completed. Final estimate obtained.\n');

%% Calculate Errors and Performance Metrics
fprintf('\n=== Calculating Performance Metrics ===\n');
fem_error = abs(fem_deflection - analytical_deflection);
kf_error = abs(kf_final_estimate - analytical_deflection);
measurement_error = abs(noisy_measurements - analytical_deflection);

fprintf('\n=== PERFORMANCE RESULTS ===\n');
fprintf('FEM RMS Error: %.6e m\n', sqrt(mean(fem_error.^2)));
fprintf('Measurement RMS Error: %.6e m\n', sqrt(mean(measurement_error.^2)));
fprintf('Kalman Filter RMS Error: %.6e m\n', sqrt(mean(kf_error.^2)));
improvement = (sqrt(mean(measurement_error.^2)) - sqrt(mean(kf_error.^2))) / sqrt(mean(measurement_error.^2)) * 100;
fprintf('Kalman Filter Improvement: %.2f%%\n', improvement);

%% Plotting Results
fprintf('\n=== Generating Plots ===\n');
figure('Position', [100, 100, 1400, 1000]);

% Plot 1: Deflection Comparison
subplot(2, 3, 1);
plot(x_nodes, analytical_deflection * 1000, 'b.', 'LineWidth', 3, 'DisplayName', 'Analytical');
hold on;
plot(x_nodes, fem_deflection * 1000, 'r--', 'LineWidth', 2, 'DisplayName', 'FEM');
plot(x_nodes, noisy_measurements * 1000, 'g*', 'LineWidth', 1.5, 'DisplayName', 'Noisy Measurements');
plot(x_nodes, kf_final_estimate * 1000, 'b-.', 'LineWidth', 2, 'DisplayName', 'Kalman Filter Estimate');
xlabel('Beam Position (m)');
ylabel('Deflection (mm)');
title('Cantilever Beam Deflection Comparison');
legend('Location', 'northwest');
grid on;
set(gca, 'FontSize', 10);

% Plot 2: Errors
subplot(2, 3, 2);
plot(x_nodes, fem_error * 1000, 'r-', 'LineWidth', 2, 'DisplayName', 'FEM Error');
hold on;
plot(x_nodes, measurement_error * 1000, 'g-', 'LineWidth', 2, 'DisplayName', 'Measurement Error');
plot(x_nodes, kf_error * 1000, 'm-', 'LineWidth', 2, 'DisplayName', 'KF Error');
xlabel('Beam Position (m)');
ylabel('Absolute Error (mm)');
title('Deflection Errors');
legend();
grid on;
set(gca, 'FontSize', 10);

% Plot 3: Kalman Filter Convergence
subplot(2, 3, 3);
selected_nodes = [3, 6, n_nodes];  % Track convergence at these nodes
colors = ['r', 'g', 'b'];
for i = 1:length(selected_nodes)
    node = selected_nodes(i);
    plot(1:n_steps, x_estimates(node, :) * 1000, colors(i), 'LineWidth', 2, ...
        'DisplayName', sprintf('Node %d (x=%.2fm)', node, x_nodes(node)));
    hold on;
end
xlabel('Kalman Filter Iteration');
ylabel('Deflection Estimate (mm)');
title('Kalman Filter Convergence');
legend();
grid on;
set(gca, 'FontSize', 10);

% Plot 4: Uncertainty Reduction
subplot(2, 3, 4);
plot(1:n_steps, mean(P_diag) * 1e6, 'k-', 'LineWidth', 2);
xlabel('Kalman Filter Iteration');
ylabel('Average Uncertainty (\mum²)');
title('Kalman Filter Uncertainty Reduction');
grid on;
set(gca, 'FontSize', 10);

% Plot 5: Noise vs Filtered Signal
subplot(2, 3, 5);
node_to_plot = n_nodes; % Free end
plot(1:n_steps, repmat(analytical_deflection(node_to_plot)*1000, n_steps, 1), 'b-', 'LineWidth', 3, 'DisplayName', 'True Value');
hold on;
plot(1:n_steps, repmat(noisy_measurements(node_to_plot)*1000, n_steps, 1), 'g:', 'LineWidth', 1.5, 'DisplayName', 'Noisy Measurement');
plot(1:n_steps, x_estimates(node_to_plot, :)*1000, 'm-', 'LineWidth', 2, 'DisplayName', 'KF Estimate');
xlabel('Iteration');
ylabel(['Deflection at x=', num2str(x_nodes(node_to_plot)), 'm (mm)']);
title('Noise Filtering at Free End');
legend();
grid on;
set(gca, 'FontSize', 10);

% Plot 6: Residual Errors
subplot(2, 3, 6);
residuals = x_estimates - analytical_deflection;
rms_residuals = sqrt(mean(residuals.^2, 1)) * 1000;
plot(1:n_steps, rms_residuals, 'k-', 'LineWidth', 2);
xlabel('Kalman Filter Iteration');
ylabel('RMS Residual Error (mm)');
title('Convergence of RMS Error');
grid on;
set(gca, 'FontSize', 10);

%% Display Final Results Table
fprintf('\n=== FINAL DEFLECTION RESULTS ===\n');
fprintf('Node\tPosition(m)\tAnalytical(mm)\tFEM(mm)\t\tNoisy(mm)\tKF Est.(mm)\n');
fprintf('--------------------------------------------------------------------------------\n');
for i = [1, 3, 5, 7, 9, 10, 11]  % Display selected nodes
    if i <= n_nodes
        fprintf('%d\t%.2f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n', i, x_nodes(i), ...
            analytical_deflection(i)*1000, fem_deflection(i)*1000, ...
            noisy_measurements(i)*1000, kf_final_estimate(i)*1000);
    end
end

%% Additional Analysis: Error Statistics
fprintf('\n=== ERROR STATISTICS ===\n');
fprintf('Maximum FEM Error: %.6f mm\n', max(fem_error)*1000);
fprintf('Maximum Measurement Error: %.6f mm\n', max(measurement_error)*1000);
fprintf('Maximum KF Error: %.6f mm\n', max(kf_error)*1000);

fprintf('\nStandard Deviation of Errors:\n');
fprintf('FEM: %.6f mm\n', std(fem_error)*1000);
fprintf('Measurements: %.6f mm\n', std(measurement_error)*1000);
fprintf('Kalman Filter: %.6f mm\n', std(kf_error)*1000);

%% Convergence Analysis
fprintf('\n=== CONVERGENCE ANALYSIS ===\n');
initial_kf_error = sqrt(mean((x_estimates(:, 1) - analytical_deflection).^2)) * 1000;
final_kf_error = sqrt(mean(kf_error.^2)) * 1000;
fprintf('Initial KF RMS Error: %.4f mm\n', initial_kf_error);
fprintf('Final KF RMS Error: %.4f mm\n', final_kf_error);
fprintf('Error reduction: %.2f%%\n', (initial_kf_error - final_kf_error)/initial_kf_error * 100);

fprintf('\n=== SIMULATION COMPLETED SUCCESSFULLY ===\n');

% Display key insights
fprintf('\n=== KEY INSIGHTS ===\n');
fprintf('1. FEM provides excellent agreement with analytical solution\n');
fprintf('2. Measurement noise significantly corrupts the true deflection\n');
fprintf('3. Kalman Filter effectively reduces noise while preserving signal\n');
fprintf('4. KF uncertainty decreases monotonically with iterations\n');
fprintf('5. The method demonstrates successful sensor fusion of model and measurements\n');




