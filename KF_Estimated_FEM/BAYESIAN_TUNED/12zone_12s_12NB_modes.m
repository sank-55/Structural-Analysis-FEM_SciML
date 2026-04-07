
%%
clear all; close all; clc;

%% 1. PHYSICAL PARAMETERS & SYSTEM SETUP
p=1000; ro=2700; yo=70*10^9; neu=0.3; h=0.001;
a=0.6; b=0.4; Nm=6; Nm_nm=12; Ndam=12; N=5000; deltat=0.0001; w_1=100*pi;
nObserv=12; D = (yo*h^3)/(12*(1-neu^2));

% %  % Damage and Sensor Locations
% % damdox = [0 a/3; a/3 2*a/3; 2*a/3 a; 0 a/3; a/3 2*a/3; 2*a/3 a];
% % damdoy = [0 b/2; 0 b/2; 0 b/2; b/2 b; b/2 b; b/2 b];
 % 12 Damage Zone Coordinates (4x3 Grid over the plate)
damdox = [0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a; ...
          0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a; ...
          0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a];

damdoy = [0 b/3; 0 b/3; 0 b/3; 0 b/3; ...
          b/3 2*b/3; b/3 2*b/3; b/3 2*b/3; b/3 2*b/3; ...
          2*b/3 b; 2*b/3 b; 2*b/3 b; 2*b/3 b];

sensox = [0.1220, 0.4780, 0.1220, 0.4780, 0.4881, 0.1119, 0.2949, 0.3051, 0.3966, 0.2034, 0.3254, 0.2136];
sensoy = [0.0718, 0.0718, 0.3282, 0.3282, 0.2359, 0.2359, 0.3282, 0.0718, 0.1538, 0.1538, 0.2359, 0.2564];
 % % sensox=[a/6  5*a/6  a/6  5*a/6   a/2 a/2];
 % % sensoy=[b/4  b/4  3*b/4  3*b/4  b/2  b/4];

%%
% SYMBOLIC EXPRESSION 
% Modal Shape and Differential Operators
syms x y
phi = [sin(pi*x/a)*sin(pi*y/b), sin(pi*x/a)*sin(2*pi*y/b), ...
       sin(2*pi*x/a)*sin(pi*y/b), sin(2*pi*x/a)*sin(2*pi*y/b), ...
       sin(1*pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(1*pi*y/b)];
dxxphi = diff(phi,x,2); dxyphi = diff(diff(phi,x),y); dyyphi = diff(phi,y,2);

phi_nm = [  sin(1*pi*x/a)*sin(1*pi*y/b), sin(1*pi*x/a)*sin(2*pi*y/b), ...
            sin(2*pi*x/a)*sin(1*pi*y/b), sin(2*pi*x/a)*sin(2*pi*y/b), ...
            sin(1*pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(1*pi*y/b), ...
            sin(2*pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(2*pi*y/b), ...
            sin(3*pi*x/a)*sin(3*pi*y/b), sin(1*pi*x/a)*sin(4*pi*y/b), ...
            sin(4*pi*x/a)*sin(1*pi*y/b), sin(2*pi*x/a)*sin(4*pi*y/b)]; 
            %sin(4*pi*x/a)*sin(2*pi*y/b), sin(3*pi*x/a)*sin(4*pi*y/b), ...
            %sin(4*pi*x/a)*sin(3*pi*y/b), sin(4*pi*x/a)*sin(4*pi*y/b)];

dxxphi_nm = diff(phi_nm,x,2); dxyphi_nm = diff(diff(phi_nm,x),y); dyyphi_nm = diff(phi_nm,y,2);

% Precompute Global Matrices
fprintf('Computing Global Matrices...\n');
Mgm = zeros(Nm); Kgk_undamaged = zeros(Nm); Fgf = zeros(Nm,1);
for i=1:Nm
    for j=1:Nm
        Mgm(i,j) = double(ro*h*int(int(phi(i)*phi(j),x,0,a),y,0,b));
        Kgk_undamaged(i,j) = double(D*int(int(dxxphi(i)*dxxphi(j)+2*dxyphi(i)*dxyphi(j)+dyyphi(i)*dyyphi(j),x,0,a),y,0,b));
    end
    Fgf(i,1) = double(int(int(p*phi(i),x,0,a),y,0,b));
end
Cgc = 0.0003*Mgm + 0.0003*Kgk_undamaged;

% Precompute Damage Influence Matrices (Kd)
Kd = zeros(Ndam, Nm, Nm);
for i=1:Ndam
    for j=1:Nm
        for k=1:Nm
            Kd(i,j,k) = double(D*int(int(dxxphi(j)*dxxphi(k)+2*dxyphi(j)*dxyphi(k)+dyyphi(j)*dyyphi(k),x,damdox(i,1),damdox(i,2)),y,damdoy(i,1),damdoy(i,2)));
        end
    end
end

% M ,K , C for nm time integration
Mgm_nm = zeros(Nm_nm); Kgk_undamaged_nm = zeros(Nm_nm); Fgf_nm = zeros(Nm_nm,1);
for i=1:Nm_nm
    for j=1:Nm_nm
        Mgm_nm(i,j) = double(ro*h*int(int(phi_nm(i)*phi_nm(j),x,0,a),y,0,b));
        Kgk_undamaged_nm(i,j) = double(D*int(int(dxxphi_nm(i)*dxxphi_nm(j)+2*dxyphi_nm(i)*dxyphi_nm(j)+dyyphi_nm(i)*dyyphi_nm(j),x,0,a),y,0,b));
    end
    Fgf_nm(i,1) = double(int(int(p*phi_nm(i),x,0,a),y,0,b));
end
Cgc_nm = 0.0003*Mgm_nm + 0.0003*Kgk_undamaged_nm;

% Precompute Damage Influence Matrices (Kd)
Kd_nm = zeros(Ndam, Nm_nm, Nm_nm);
for i=1:Ndam
    for j=1:Nm_nm
        for k=1:Nm_nm
            Kd_nm(i,j,k) = double(D*int(int(dxxphi_nm(j)*dxxphi_nm(k)+2*dxyphi_nm(j)*dxyphi_nm(k)+dyyphi_nm(j)*dyyphi_nm(k),x,damdox(i,1),damdox(i,2)),y,damdoy(i,1),damdoy(i,2)));
        end
    end
end

% Sensor Observation Matrix (Phi)
Phi_mat = zeros(nObserv, Nm);
for i=1:nObserv
    x1=sensox(i); y1=sensoy(i);
    Phi_mat(i,:) = [sin(pi*x1/a)*sin(pi*y1/b), sin(1*pi*x1/a)*sin(2*pi*y1/b), ...
                    sin(2*pi*x1/a)*sin(1*pi*y1/b), sin(2*pi*x1/a)*sin(2*pi*y1/b), ...
                    sin(pi*x1/a)*sin(3*pi*y1/b), sin(3*pi*x1/a)*sin(pi*y1/b)];
end

Phi_mat_nm = zeros(nObserv, Nm_nm);
for i=1:nObserv
    x1=sensox(i); y1=sensoy(i);
    Phi_mat_nm(i,:) = [ sin(1*pi*x1/a)*sin(1*pi*y1/b), sin(1*pi*x1/a)*sin(2*pi*y1/b), ...
                        sin(2*pi*x1/a)*sin(1*pi*y1/b), sin(2*pi*x1/a)*sin(2*pi*y1/b), ...
                        sin(1*pi*x1/a)*sin(3*pi*y1/b), sin(3*pi*x1/a)*sin(1*pi*y1/b), ...
                        sin(2*pi*x1/a)*sin(3*pi*y1/b), sin(3*pi*x1/a)*sin(2*pi*y1/b), ...
                        sin(3*pi*x1/a)*sin(3*pi*y1/b), sin(1*pi*x1/a)*sin(4*pi*y1/b), ...
                        sin(4*pi*x1/a)*sin(1*pi*y1/b), sin(2*pi*x1/a)*sin(4*pi*y1/b)]; 
                        %sin(4*pi*x1/a)*sin(2*pi*y1/b), sin(3*pi*x1/a)*sin(4*pi*y1/b), ...
                        %sin(4*pi*x1/a)*sin(3*pi*y1/b), sin(4*pi*x1/a)*sin(4*pi*y1/b)];
end



%% SAVIG THE SYM EXP
save('PlateModelData_12NM.mat', 'Mgm', 'Kgk_undamaged', 'Kd', 'Phi_mat','Mgm_nm', 'Kgk_undamaged_nm', 'Kd_nm', 'Phi_mat_nm', 'Fgf', 'Fgf_nm');

%% LOADING THE DATA 
load('PlateModelData_12NM.mat');
%% 2. GENERATE TRUTH DATA (Damaged Simulation)
% damf_true = zeros(Ndam, N); X_DAM =linspace(0, 0.8, N-199);
% damf_true(2, 1200:N) = 0.6 * sqrt( (0:N-1200)/(N-1200) );
% damf_true(2, 1200:N) = linspace(0, 0.6, N-1199);
% damf_true(5, 30:N) = linspace(0, 0.2, N-29);

damf_true=zeros(Ndam,N);
% 6 dam zones 
% % 
%   for i=200:N
%       damf_true(1,i)=(i-200)^.7*0.3/(N-200)^.7;
%   end
% % % 
%   for i=1200:N
%       damf_true(2,i)=(i-1200)*0.6/(N-1200);
%   end
% 
% %   for i=600:N
% %      damf_true(3,i)=(i-600)*0.12/(N-600);
% %  end
% % 
% % for i=10:N
% %     damf_true(4,i)=((i-10)^0.8)*0.3/(N-10);
% % end
% 
% for i=30:N
%     damf_true(5,i)=(i-30)*0.2/(N-30);
% end
% 
% 
% for i=300:N
%     damf_true(6,i)=(i-300)*0.2/(N-300);
% end

% for 12 zones 
for i = 200:N,  damf_true(1,i) = (i-200)^.7 * 0.3 / (N-200)^.7; end
for i = 1200:N, damf_true(2,i) = (i-1200) * 0.6 / (N-1200); end
for i = 600:N,  damf_true(3,i) = (i-600) * 0.3 / (N-600); end
for i = 30:N,   damf_true(5,i) = (i-30) * 0.2 / (N-30); end
for i = 800:N,  damf_true(12,i) = (i-800) * 0.1 / (N-800); end
for i = 1500:N, damf_true(11,i) = (i-1500) * 0.25 / (N-1500); end


% Run Newmark-Beta to get "Measured" Sensors

gamma = 0.5; 
beta = 0.25;

% Precompute Newmark Constants
a0 = 1/(beta*deltat^2); a1 = gamma/(beta*deltat);
a2 = 1/(beta*deltat);   a3 = 1/(2*beta) - 1;
a4 = gamma/beta - 1;    a5 = deltat/2 * (gamma/beta - 2);

% Initial Conditions
x_t = zeros(Nm_nm, 1); 
v_t = zeros(Nm_nm, 1); 
% Solve for initial acceleration: a = M \ (F0 - C*v0 - K*x0)
a_t = Mgm_nm \ (Fgf_nm*sin(w_1*0) - Cgc_nm*v_t - Kgk_undamaged_nm*x_t); 
zst_measured = zeros(2*nObserv, N);
zacc_measured = zeros(nObserv, N);

for i = 2:N
    t = (i-1)*deltat;
    
    % 1. Update current Damage/Stiffness
    K_curr = Kgk_undamaged_nm;
    for d = 1:Ndam
        K_curr = K_curr - damf_true(d,i) * squeeze(Kd_nm(d,:,:)); 
    end
    
    % 2. Effective Stiffness Matrix
    K_hat = K_curr + a1*Cgc_nm + a0*Mgm_nm;
    
    % 3. Effective Force Vector (Predictor)
    F_ext = Fgf_nm * sin(w_1*t);
    F_hat = F_ext + Mgm_nm*(a0*x_t + a2*v_t + a3*a_t) + Cgc_nm*(a1*x_t + a4*v_t + a5*a_t);
    
    % 4. Solve for New Displacement
    x_new = K_hat \ F_hat;
    
    % 5. Derived Velocity and Acceleration (Corrector)
    a_new = a0*(x_new - x_t) - a2*v_t - a3*a_t;
    v_new = v_t + (1-gamma)*deltat*a_t + gamma*deltat*a_new;
    
    % 6. Store Observation and update states
    zst_measured(:,i) = [Phi_mat_nm *( x_new + x_new*0.02*randn); Phi_mat_nm *( v_new + v_new*0.02*randn)];
    zacc_measured(:,i) = Phi_mat_nm *( a_new + a_new*0.03*randn);
    
    x_t = x_new;
    v_t = v_new;
    a_t = a_new;
end





% storing variables for the DEKF







%% 3. BAYESIAN OPTIMIZATION

optimVars_ekf = [
    optimizableVariable('p_state_exp', [-19, -12], 'Type', 'real')
    optimizableVariable('p_param_exp', [-19, -17], 'Type', 'real')
    optimizableVariable('q_state_exp', [-15, -12], 'Type', 'real')
    optimizableVariable('q_param_exp', [-10, -7], 'Type', 'real')
    %optimizableVariable('p_exp', [-10, -6], 'Type', 'real')
    %optimizableVariable('q_exp', [-20, -12], 'Type', 'real')
    optimizableVariable('r_exp', [-4, 2], 'Type', 'real')
];

ObjFcn_ekf = @(p) run_ekf_internal(p, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                               nObserv, Nm, Ndam, deltat, N, w_1, zacc_measured, damf_true);



optimVars_dkf2 = [
    optimizableVariable('p_state_dkf2', [-24, -20], 'Type', 'real')
    optimizableVariable('p_param_dkf2', [-22, -18], 'Type', 'real')
    optimizableVariable('q_state_dkf2', [2, 6], 'Type', 'real')
    optimizableVariable('q_param_dkf2', [-8, -4], 'Type', 'real')
    optimizableVariable('r_state_dkf2', [-1, 3], 'Type', 'real')
    optimizableVariable('r_param_dkf2', [1, 4], 'Type', 'real')
    optimizableVariable('rval', [0.2, 0.9], 'Type', 'real')

    ];
 ObjFcn_dkf2 = @(p) run_dkf2_internal(p, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);


optimVars_dkf1 = [
    optimizableVariable('p_state_dkf1', [-22, -20], 'Type', 'real')
    optimizableVariable('p_param_dkf1', [-22, -18], 'Type', 'real')
    optimizableVariable('q_state_dkf1', [8, 12], 'Type', 'real')
    optimizableVariable('q_param_dkf1', [-12, -8], 'Type', 'real')
    optimizableVariable('r_state_dkf1', [1, 4], 'Type', 'real')
    optimizableVariable('r_param_dkf1', [-2, 2], 'Type', 'real')
    
    ];
 ObjFcn_dkf1 = @(p) run_dkf1_internal(p, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);


fprintf('\nStarting Bayesian Optimization ekf ...\n');
results_ekf = bayesopt(ObjFcn_ekf, optimVars_ekf, 'MaxObjectiveEvaluations', 200, 'PlotFcn', {});

fprintf('\nStarting Bayesian Optimization dkf1 ...\n');
results_dkf1 = bayesopt(ObjFcn_dkf1, optimVars_dkf1, 'MaxObjectiveEvaluations', 200, 'PlotFcn', {});

fprintf('\nStarting Bayesian Optimization dkf2 ...\n');
results_dkf2 = bayesopt(ObjFcn_dkf2, optimVars_dkf2, 'MaxObjectiveEvaluations', 200, 'PlotFcn', {});

%% 4. FINAL SIMULATION WITH BEST PARAMETERS
bestP = bestPoint(results_ekf);
[final_rmse, fd_est ] = run_ekf_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1, zacc_measured, damf_true);

bestP = bestPoint(results_dkf1);
% for scheme 2
[final_rmse_dkf1, fd_dkf1_optimized ] = run_dkf1_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);

bestP = bestPoint(results_dkf2);
% for scheme 2
[final_rmse_dkf2, fd_dkf2_optimized ] = run_dkf2_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);

%%
figure(1);
for i=1:Ndam
    subplot(3,4,i);
    plot(damf_true(i,:), 'k', 'LineWidth', 1.5); hold on;
    %plot(fd_dkf2_optimized(i,:), 'b--', 'LineWidth', 1);
    plot(fd_est(i,:), 'g--', 'LineWidth', 1);
    title(['Zone ', num2str(i)]); grid on;
    if i==1; legend('True','EKF Tuned'); end
end
%%
figure(17);
clf; 
set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1, 1, 12, 9]); 

marker_interval = 450; 
smooth_span = 100;
smooth_span_ekf = 100;
for i = 1:Ndam
    % --- 1. Calculate RMSE ---
    rmse1 = sqrt(mean((damf_true(i,:) - fd_dkf1_optimized(i,:)).^2));
    rmse2 = sqrt(mean((damf_true(i,:) - fd_dkf2_optimized(i,:)).^2));
    rmse3 = sqrt(mean((damf_true(i,:) - fd_est(i,:)).^2));
    
    subplot(3,4,i);
    hold on;
    
    % --- 2. Plotting with "Haloed" Markers ---
    h_true = plot(damf_true(i,:), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 2.5);
    
    % DKF1
    smooth_dk1 = movmean(fd_dkf1_optimized(i,:), smooth_span);
    h_dk1 = plot(smooth_dk1, 'k-', 'LineWidth', 1.1, ...
        'Marker', 'o', 'MarkerIndices', 1 : marker_interval : length(smooth_dk1), ...
        'MarkerSize', 4.5, 'MarkerFaceColor', 'w');
    
    % DKF2
    smooth_dk2 = movmean(fd_dkf2_optimized(i,:), smooth_span);
    h_dk2 = plot(smooth_dk2, 'k-', 'LineWidth', 1.1, ...
        'Marker', 's', 'MarkerIndices', floor(marker_interval/3) : marker_interval : length(smooth_dk2), ...
        'MarkerSize', 4.5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
    
    % EKF
    smooth_ekf = movmean(fd_est(i,:), smooth_span_ekf);
    h_ekf = plot(smooth_ekf, 'k:.', 'LineWidth', 1.1, ...
        'Marker', '+', 'MarkerIndices', floor(2*marker_interval/3) : marker_interval : length(smooth_ekf), ...
        'MarkerSize', 6);

    % --- 3. Refined RMSE Text Placement ---
    % Placed at top-left (0.05, 0.9) of each subplot
    text_str = {sprintf('DK1: %.4f', rmse1), ...
                sprintf('DK2: %.4f', rmse2), ...
                sprintf('EKF: %.4f', rmse3)};
    text(0.05, 0.90, text_str, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 7.5, ...
        'FontName', 'Helvetica', 'BackgroundColor', [1 1 1 0.7]);

    % --- 4. Subplot Formatting ---
    title(['Zone ', num2str(i)], 'FontSize', 10, 'FontWeight', 'bold');
    grid on; box on;
    set(gca, 'GridAlpha', 0.1, 'FontSize', 9);
    if mod(i-1, 4) == 0, ylabel('Damage Factor'); end
    if i > 8, xlabel('Time Step'); end
    
    if i == 1, plot_handles = [h_true, h_dk1, h_dk2, h_ekf]; end
end

% --- 5. Corrected Legend Placement ---
% [left, bottom, width, height] in figure-normalized units
% Values around 0.5 center the legend horizontally
leg_labels = {'True', 'DKF1[RMSE:0.0214]', 'DKF2[RMSE:0.0203]', ...
              'EKF[RMSE:0.0374]  \bf[Sensor:12]\rm'};
L = legend(plot_handles, leg_labels, 'Orientation', 'horizontal', ...
    'Interpreter', 'tex', 'FontSize', 11);
set(L, 'Position', [0.15, 0.015, 0.7, 0.03], 'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);

%% --- HELPER: EKF FUNCTION ---
function [score, fd_history] = run_ekf_internal(p, M, K0, C, Kd_all, Phi, Fgf, nObs, Nm, Ndam, dt, N, w1, z, truth)
    % Unpack
    P_state = 10^p.p_state_exp;
    P_param = 10^p.p_param_exp;
    Q_state = 10^p.q_state_exp;
    Q_param = 10^p.q_param_exp;
    %P_val = 10^(p.p_exp);
    %Q_val = 10^(p.q_exp);
    R_val = 10^p.r_exp;
    state_ekf_store = zeros(Nm,N);
    
    % Init
    xx = zeros(2*Nm + Ndam, 1);
    %P = blkdiag(1e-6*eye(2*Nm), 1e-6*eye(Ndam));
    P = 1*blkdiag(P_state*eye(2*Nm), P_param*eye(Ndam));
    Q = 1*blkdiag(Q_state*eye(2*Nm), Q_param*eye(Ndam));
    R = R_val * eye(nObs);
    fd_history = zeros(Ndam, N);
    Obsv_ekf= zeros(2*Nm+Ndam,1);





    for i=2:N
        u = xx(1:Nm); v = xx(Nm+1:2*Nm); fd = xx(2*Nm+1:end);
        
        % Jacobian Prep
        K_eff = K0; S = zeros(Nm, Ndam);
        for d=1:Ndam
            Ki = squeeze(Kd_all(d,:,:));
            K_eff = K_eff - fd(d)*Ki;
            S(:,d) = Ki * u;
        end
        
        % Continuous Jacobian A
        A = [zeros(Nm), eye(Nm), zeros(Nm,Ndam);
             -M\K_eff, -M\C, M\S;
             zeros(Ndam, 2*Nm+Ndam)];
         
        % 2nd Order Taylor Discretization: F = I + Adt + 0.5(Adt)^2
        At = A*dt;
        F_jac = eye(2*Nm+Ndam) + At + 0.5*(At*At);
        
        % Prediction
        acc = M \ (Fgf*sin(w1*(i-1)*dt) - K_eff*u - C*v);
        xx(1:Nm) = u + v*dt;
        xx(Nm+1:2*Nm) = v + acc*dt;
        P = F_jac * P * F_jac' + Q;
        
        % Update
        H = [-Phi*(M\K_eff), -Phi*(M\C), Phi*(M\S)];
        innov = z(:,i) - Phi*acc;
        K_gain = (P*H') / (H*P*H' + R);
        xx = xx + K_gain * innov;
        P = (eye(length(xx)) - K_gain*H) * P;
        
        state_ekf_store(:,i)=xx(1:Nm);
        ob_ekf = obsv(F_jac,H);


         for k=1: 2*Nm+Ndam

            Obsv_ekf(k) = Obsv_ekf(k) + norm(ob_ekf(:, k))^2;
        
        end



        fd_history(:,i) = xx(2*Nm+1:end);
    end
    figure(10);
    bar(Obsv_ekf(2*Nm+1:2*Nm+Ndam));%/max(Obsv_fd_s1));
    set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6','Zone 7', 'Zone 8', 'Zone 9', 'Zone 10', 'Zone 11', 'Zone 12'});
    title('Observability(Fd) over Time for EKF');
    
    score = sqrt(mean((truth(:) - fd_history(:)).^2));
end




%% DEKF1 loop 
function[final_rmse_dkf1, fd_dkf1_opt ] = run_dkf1_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf)

X=zeros(nObserv,Nm);


H=[Phi_mat  X ; X  Phi_mat];



fac=zeros(Ndam,1);
facc=fac;
fd_ex = zeros(Ndam,1);
fd_est_ex = zeros(Ndam,N);
fd_dkf2 = zeros(Ndam,1);
fd_est_dkf2 = zeros(Ndam,N);


xxp=zeros(2*Nm,1);
xxp2=zeros(2*Nm,1);
xxp_exfd=zeros(2*Nm+Ndam,1);  % for the extended kalman filter 
vvp=zeros(Nm,1);
aap=zeros(Nm,1);
accp=zeros(Nm,1);
   acp=zeros(Nm,1);
   acp_2 = zeros(Nm,1);
   xp=zeros(Nm,1);
   vp=zeros(Nm,1);
   xrip=zeros(Nm,1);
   vrip=zeros(Nm,1);
   accpri=zeros(Nm,1);
xx_true_store = zeros(Nm,N);   
xx_store_dkf1 = zeros(Nm,N);
xx_store_dkf2 = zeros(Nm,N);
xx_store_ekf = zeros(Nm,N);
%% Creating some obs matrix for the different schemes and the states ( sc1 & sc2 & extened) 
Obsv_state_s1 = zeros(2*Nm,1);
Obsv_fd_s2 = zeros(Ndam,1);
Obsv_fd_s1 = zeros(Ndam,1);
Obsv_state_s2 = zeros(2*Nm,1);
Obsv_ekf = zeros(2*Nm + Ndam, 1);


% params for dekf
Pdkf2 = eye(2*Nm)*10^(bestP.p_state_dkf1);
%Pdkf2 = eye(2*Nm)*10^(-2);
Qdkf2 = eye(2*Nm)*10^(bestP.q_state_dkf1);
Pf2 = blkdiag(0.1288,0.0468,0.0494,0.1424,0.1454, 0.0150,0.0151,0.1247,0.1187,0.0621,0.0638,0.1456)*10^(bestP.p_param_dkf1);
%Pf2 = 1e-5*blkdiag(0.1077 , 0.0544 ,0.1292, 0.1234 ,  0.0522,  0.1149);
Qf2 = eye(Ndam)*10^(bestP.q_param_dkf1);
Rdkf2 = eye(2*nObserv)*10^(bestP.r_state_dkf1);
Rf2 = eye(nObserv)*10^(bestP.r_param_dkf1);
%r = bestP.rval;

%  P=[ones(Nm,1)*10^(-2);ones(Nm,1)*10^(-2)];   % 10^-10  for in 4 sensors
%  P=diag(P);
%  Pdkf2=P;
% 
%  Q=eye(2*Nm)*10^-16;  %  % scheme 1 -> eye(dim)/10^-10 for 6 sensor (unchange modes)
%  Qdkf2=eye(2*Nm)*10^-16; % scheme 2 -> eye(dim)/10^-8 for 6 sensor (unchange modes)
% 
% Pf = 1e-5*blkdiag( 0.2135  , 0.1065 , 0.2160, 0.2371, 0.1551 ,  0.2345);% for 6 sensor scheme 1;
% Pf2 = 1e-5*blkdiag(0.7579 ,  0.7070 ,0.6272, 0.6903 ,  0.7201,  0.6761); % for 6 sensor scheme 2
% Qf =  1e-6*Pf; %1e-12*blkdiag( 0.4793  , 0.3464 ,0.5248, 0.4994 ,0.3564, 0.5304);
% Qf2 = 1e-6*Pf2; %    % 1e1*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);
% 
% R=1*eye(2*nObserv)*10^(0.34927);  % for state ( disp & vel) of Scheme - 1 for 6 sen ->1*eye(2*nObserv)*10^(-4) works better ( no noise)  
% Rf=1*eye(nObserv)*10^(0.34927);
% Rf2=1*eye(nObserv)*10^(0.34927);
% %Rf = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);
% 
% Rdkf2=1*eye(2*nObserv)*10^(0.34927); % for state ( disp & vel) of Scheme - 2
%Rf2 = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);

for i =2:N
    t = (i-1)*deltat;
    F=Fgf*sin(w_1*t);    
    F1=F;
    Mgminv=inv(Mgm);
    Af=eye(Ndam);

    Kgks2 = Kgk_undamaged;
     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgks2=Kgks2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
     end

      Kdx2=zeros(Nm,Ndam);

    for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        xi(:,j)=tem*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) );
        Kdx2(:,j)=xi(:,j);
    end
     
    Hdkf2 = Phi_mat*(inv(Mgm)*(Kdx2)) ;
   
    a_dkf2 = Hdkf2*fd_dkf2 + Phi_mat*(inv(Mgm)*(F1- Kgk_undamaged*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm)) - Cgc*(xxp2(Nm+1:2*Nm) + deltat*acp_2) ));
     

    Pnf2=Af*Pf2*transpose(Af)+Qf2;
    K_gain_dkf2=Pnf2*transpose(Hdkf2)*inv((Hdkf2*Pnf2*transpose(Hdkf2))+Rf2);
    
    %Af=eye(1);
    reba=100;
    fd_dkf2=fd_dkf2+K_gain_dkf2*(zacc_measured(:,i) - a_dkf2);



      Kgdkf2 = Kgk_undamaged;
     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgdkf2=Kgdkf2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
     end


        Adkf2=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgdkf2  -inv(Mgm)*Cgc];
        Bdkf2=[zeros(Nm,Nm); inv(Mgm)];

    

    Ax2 = eye(2*Nm) + deltat*Adkf2 + 0.5*deltat^2*(Adkf2*Adkf2);
    Bx2 = deltat*Bdkf2 ;

    % % How this is defined
    % Ax2=deltat*Adkf2+eye(2*Nm);
    % Bx2=deltat*Bdkf2;
    % Ax2=zeros(2*Nm,2*Nm);
    % Bx2_=zeros(2*Nm,2*Nm);
    % Bx2=zeros(2*Nm,Nm);
    % for jj=0:100
    %     Ax2=Ax2+(Adkf2*deltat)^jj/factorial(jj);
    %     Bx2_=Bx2_+(Adkf2*deltat)^jj/factorial(jj+1);
    % end

    %Bx2=Bx2_*Bdkf2*deltat;
        % For Scheme 2 
    xx2=Ax2*xxp2+Bx2*F;
    Pndkf2=Ax2*Pdkf2*transpose(Ax2)+Qdkf2;    
    Kx2=Pndkf2*transpose(H)*inv((H*Pndkf2*transpose(H))+Rdkf2);
    xx2=xx2+Kx2*(zst_measured(:,i)-H*xx2);   
    xx_store_dkf2(:,i) = xx2(1:Nm);
    fd_est_dkf2(:,i) = fd_dkf2;

    O_s1_st = obsv(Adkf2,H);
    O_s1_fd = obsv(Af,Hdkf2);

    acp_2 = (xx2(Nm+1:2*Nm)-xxp2(Nm+1:2*Nm))/deltat;
    xxp2 = xx2;
    Pdkf2=(eye(2*Nm)-Kx2*H)*Pndkf2;
    Pf2=(eye(Ndam)-K_gain_dkf2*Hdkf2)*Pnf2;

    

    % updating the obsv terms
    for z=1:Ndam
        Obsv_fd_s1(z) = Obsv_fd_s1(z) + norm(O_s1_fd(:, z))^2;

   
    end

    for k=1: 2*Nm

        Obsv_state_s1(k) = Obsv_state_s1(k) + norm(O_s1_st(:, k))^2;
        
    end

end




%% Plotting the observation matrix 
figure(6)
bar(Obsv_fd_s1);%/max(Obsv_fd_s1));
set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6','Zone 7', 'Zone 8', 'Zone 9', 'Zone 10', 'Zone 11', 'Zone 12'});
title('Observability(Fd) over Time for scheme1');



figure(8)
bar(Obsv_state_s1);%/max(Obsv_state_s1));
set(gca, 'XTickLabel', {'m 1', 'm 2', 'm 3', 'm 4', 'm 5', 'm 6', 'm 7','m 8','m9','m10','m11','m12'});
title(' Observability Strength (modes) over Time for scheme1 ');





P1_s1 = 1./Obsv_state_s1;
P2_s1 = 1./Obsv_fd_s1;

for i=1:Ndam
    P2_s1(i);
end

%P2_s1 = 1./Obsv_fd_s1;

%P1_s1 = 1./Obsv_state_s1;
fd_dkf1_opt = fd_est_dkf2;
final_rmse_dkf1 = sqrt(mean((damf(:) - fd_est_dkf2(:)).^2));
end



%% DEKF2 FUNCTION 
function[final_rmse_dkf2, fd_dkf2_opt ] = run_dkf2_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf)

X=zeros(nObserv,Nm);


H=[Phi_mat  X ; X  Phi_mat];



fac=zeros(Ndam,1);
facc=fac;
fd_ex = zeros(Ndam,1);
fd_est_ex = zeros(Ndam,N);
fd_dkf2 = zeros(Ndam,1);
fd_est_dkf2 = zeros(Ndam,N);


xxp=zeros(2*Nm,1);
xxp2=zeros(2*Nm,1);
xxp_exfd=zeros(2*Nm+Ndam,1);  % for the extended kalman filter 
vvp=zeros(Nm,1);
aap=zeros(Nm,1);
accp=zeros(Nm,1);
   acp=zeros(Nm,1);
   acp_2 = zeros(Nm,1);
   xp=zeros(Nm,1);
   vp=zeros(Nm,1);
   xrip=zeros(Nm,1);
   vrip=zeros(Nm,1);
   accpri=zeros(Nm,1);
xx_true_store = zeros(Nm,N);   
xx_store_dkf = zeros(Nm,N);
xx_store_dkf2 = zeros(Nm,N);
xx_store_ekf = zeros(Nm,N);
%% Creating some obs matrix for the different schemes and the states ( sc1 & sc2 & extened) 
Obsv_state_s1 = zeros(2*Nm,1);
Obsv_fd_s2 = zeros(Ndam,1);
Obsv_fd_s1 = zeros(Ndam,1);
Obsv_state_s2 = zeros(2*Nm,1);
Obsv_ekf = zeros(2*Nm + Ndam, 1);


% params for dekf
Pdkf2 = eye(2*Nm)*10^(bestP.p_state_dkf2);
%Pdkf2 = eye(2*Nm)*10^(-2);
Qdkf2 = eye(2*Nm)*10^(bestP.q_state_dkf2);
Pf2 = blkdiag(0.1301,0.04768,0.0502,0.1438,0.1484, 0.0154,0.0155,0.1277,0.1197,0.0631,0.0649,0.1470)*10^(bestP.p_param_dkf2);
%Pf2 = 1e-5*blkdiag(0.1077 , 0.0544 ,0.1292, 0.1234 ,  0.0522,  0.1149);
Qf2 = eye(Ndam)*10^(bestP.q_param_dkf2);
Rdkf2 = eye(2*nObserv)*10^(bestP.r_state_dkf2);
Rf2 = eye(nObserv)*10^(bestP.r_param_dkf2);
r = bestP.rval;

%  P=[ones(Nm,1)*10^(-2);ones(Nm,1)*10^(-2)];   % 10^-10  for in 4 sensors
%  P=diag(P);
%  Pdkf2=P;
% 
%  Q=eye(2*Nm)*10^-16;  %  % scheme 1 -> eye(dim)/10^-10 for 6 sensor (unchange modes)
%  Qdkf2=eye(2*Nm)*10^-16; % scheme 2 -> eye(dim)/10^-8 for 6 sensor (unchange modes)
% 
% Pf = 1e-5*blkdiag( 0.2135  , 0.1065 , 0.2160, 0.2371, 0.1551 ,  0.2345);% for 6 sensor scheme 1;
% Pf2 = 1e-5*blkdiag(0.7579 ,  0.7070 ,0.6272, 0.6903 ,  0.7201,  0.6761); % for 6 sensor scheme 2
% Qf =  1e-6*Pf; %1e-12*blkdiag( 0.4793  , 0.3464 ,0.5248, 0.4994 ,0.3564, 0.5304);
% Qf2 = 1e-6*Pf2; %    % 1e1*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);
% 
% R=1*eye(2*nObserv)*10^(0.34927);  % for state ( disp & vel) of Scheme - 1 for 6 sen ->1*eye(2*nObserv)*10^(-4) works better ( no noise)  
% Rf=1*eye(nObserv)*10^(0.34927);
% Rf2=1*eye(nObserv)*10^(0.34927);
% %Rf = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);
% 
% Rdkf2=1*eye(2*nObserv)*10^(0.34927); % for state ( disp & vel) of Scheme - 2
%Rf2 = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);

for i =2:N
    t = (i-1)*deltat;
    F=Fgf*sin(w_1*t);    
    F1=F;
    Mgminv=inv(Mgm);
    Af=eye(Ndam);

    Kgks2 = Kgk_undamaged;
     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgks2=Kgks2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
     end

      Kdx2=zeros(Nm,Ndam);

    for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        xi(:,j)=tem*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*r*(1-r)*acp_2 );
        Kdx2(:,j)=xi(:,j);
    end
     
    Hdkf2 = Phi_mat*(inv(Mgm + deltat*(1-r)*Cgc + deltat^2*(1-r)^2*(Kgks2))*(Kdx2)) ;
   
    a_dkf2 = Hdkf2*fd_dkf2 + Phi_mat*(inv(Mgm + deltat*(1-r)*Cgc + deltat^2*(1-r)^2*(Kgks2))*(F1- Kgk_undamaged*(xxp2(1:Nm) + deltat*r*xxp2(Nm+1:2*Nm) + deltat^2*r*(1-r)*acp_2) - Cgc*(xxp2(Nm+1:2*Nm) + deltat*r*acp_2) ));
     

    Pnf2=Af*Pf2*transpose(Af)+Qf2;
    K_gain_dkf2=Pnf2*transpose(Hdkf2)*inv((Hdkf2*Pnf2*transpose(Hdkf2))+Rf2);
    
    %Af=eye(1);
    reba=100;
    fd_dkf2=fd_dkf2+K_gain_dkf2*(zacc_measured(:,i) - a_dkf2);



      Kgdkf2 = Kgk_undamaged;
     for j=1:Ndam
        for k=1:Nm
            for m=1:Nm
                tem(k,m)=Kd(j,k,m);
            end
        end
        Kgdkf2=Kgdkf2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
     end


        Adkf2=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgdkf2  -inv(Mgm)*Cgc];
        Bdkf2=[zeros(Nm,Nm); inv(Mgm)];

    

    Ax2 = eye(2*Nm) + deltat*Adkf2 + 0.5*deltat^2*(Adkf2*Adkf2);
    Bx2 = deltat*Bdkf2 ;

    % % How this is defined
    % Ax2=deltat*Adkf2+eye(2*Nm);
    % Bx2=deltat*Bdkf2;
    % Ax2=zeros(2*Nm,2*Nm);
    % Bx2_=zeros(2*Nm,2*Nm);
    % Bx2=zeros(2*Nm,Nm);
    % for jj=0:100
    %     Ax2=Ax2+(Adkf2*deltat)^jj/factorial(jj);
    %     Bx2_=Bx2_+(Adkf2*deltat)^jj/factorial(jj+1);
    % end

    %Bx2=Bx2_*Bdkf2*deltat;
        % For Scheme 2 
    xx2=Ax2*xxp2+Bx2*F;
    Pndkf2=Ax2*Pdkf2*transpose(Ax2)+Qdkf2;    
    Kx2=Pndkf2*transpose(H)*inv((H*Pndkf2*transpose(H))+Rdkf2);
    xx2=xx2+Kx2*(zst_measured(:,i)-H*xx2);   
    xx_store_dkf2(:,i) = xx2(1:Nm);
    fd_est_dkf2(:,i) = fd_dkf2;

   
    O_s2_st = obsv(Adkf2,H);
    O_s2_fd = obsv(Af,Hdkf2);

    acp_2 = (xx2(Nm+1:2*Nm)-xxp2(Nm+1:2*Nm))/deltat;
    xxp2 = xx2;
    Pdkf2=(eye(2*Nm)-Kx2*H)*Pndkf2;
    Pf2=(eye(Ndam)-K_gain_dkf2*Hdkf2)*Pnf2;
    % updating the obsv terms
    for z=1:Ndam
        %Obsv_fd_s1(z) = Obsv_fd_s1(z) + norm(O_s1_fd(:, z))^2;

        Obsv_fd_s2(z) = Obsv_fd_s2(z) + norm(O_s2_fd(:, z))^2;
    end

    for k=1: 2*Nm

        %Obsv_state_s1(k) = Obsv_state_s1(k) + norm(O_s1_st(:, k))^2;
        Obsv_state_s2(k) = Obsv_state_s2(k) + norm(O_s2_st(:, k))^2;
    end

end




%% Plotting the observation matrix 
% figure(6)
% bar(Obsv_fd_s1);%/max(Obsv_fd_s1));
% set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6'});
% title('Observability(Fd) over Time for scheme1');

figure(7)
bar(Obsv_fd_s2);%/max(Obsv_fd_s2));
set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6','Zone 7', 'Zone 8', 'Zone 9', 'Zone 10', 'Zone 11', 'Zone 12'});
title('Observability(Fd) over Time for scheme2 ');

% figure(8)
% bar(Obsv_state_s1);%/max(Obsv_state_s1));
% set(gca, 'XTickLabel', {'m 1', 'm 2', 'm 3', 'm 4', 'm 5', 'm 6', 'm 7','m 8','m9','m10','m11','m12'});
% title(' Observability Strength (modes) over Time for scheme1 ');

figure(9)
bar(Obsv_state_s2);%/max(Obsv_state_s2));
set(gca, 'XTickLabel', {'m 1', 'm 2', 'm 3', 'm 4', 'm 5', 'm 6', 'm 7','m 8','m9','m10','m11','m12'});
title(' Observability Strength (modes) over Time for scheme2 ');



P1_s2 = 1./Obsv_state_s2;
P2_s2 = 1./Obsv_fd_s2;

for i=1:Ndam
    P2_s2(i);
end

%P2_s1 = 1./Obsv_fd_s1;

%P1_s1 = 1./Obsv_state_s1;
fd_dkf2_opt = fd_est_dkf2;
final_rmse_dkf2 = sqrt(mean((damf(:) - fd_est_dkf2(:)).^2));
end
