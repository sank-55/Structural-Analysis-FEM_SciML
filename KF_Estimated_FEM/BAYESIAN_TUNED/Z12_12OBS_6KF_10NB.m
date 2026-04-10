% 
% %%
% clear all; close all; clc;
% % wORKING BUT HAVE TO CHECK FOR THE DKF2 IS WORKING BETTER OR NOT 
% 
% %% 1. PHYSICAL PARAMETERS & SYSTEM SETUP
% p=1000; ro=2700; yo=70*10^9; neu=0.3; h=0.001;
% a=0.6; b=0.4; Nm=8; Nm_nm=16; Ndam=12; N=5000; deltat=0.0001; w_1=100*pi;
% nObserv=10; D = (yo*h^3)/(12*(1-neu^2));
% 
% % %  % Damage and Sensor Locations
% % % damdox = [0 a/3; a/3 2*a/3; 2*a/3 a; 0 a/3; a/3 2*a/3; 2*a/3 a];
% % % damdoy = [0 b/2; 0 b/2; 0 b/2; b/2 b; b/2 b; b/2 b];
%  % 12 Damage Zone Coordinates (4x3 Grid over the plate)
% damdox = [0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a; ...
%           0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a; ...
%           0 a/4; a/4 a/2; a/2 3*a/4; 3*a/4 a];
% 
% damdoy = [0 b/3; 0 b/3; 0 b/3; 0 b/3; ...
%           b/3 2*b/3; b/3 2*b/3; b/3 2*b/3; b/3 2*b/3; ...
%           2*b/3 b; 2*b/3 b; 2*b/3 b; 2*b/3 b];
% 
% sensox = [0.1220, 0.4780, 0.1220, 0.4780, 0.4881, 0.1119, 0.2949, 0.3051, 0.3966, 0.2034, 0.3254, 0.2136];
% sensoy = [0.0718, 0.0718, 0.3282, 0.3282, 0.2359, 0.2359, 0.3282, 0.0718, 0.1538, 0.1538, 0.2359, 0.2564];
%  % % sensox=[a/6  5*a/6  a/6  5*a/6   a/2 a/2];
%  % % sensoy=[b/4  b/4  3*b/4  3*b/4  b/2  b/4];
% 
% %%
% % SYMBOLIC EXPRESSION 
% % Modal Shape and Differential Operators
% syms x y
% phi = [sin(pi*x/a)*sin(pi*y/b), sin(pi*x/a)*sin(3*pi*y/b), ...
%        sin(3*pi*x/a)*sin(pi*y/b), sin(3*pi*x/a)*sin(3*pi*y/b), ...
%        sin(1*pi*x/a)*sin(5*pi*y/b), sin(5*pi*x/a)*sin(1*pi*y/b),...
%        sin(1*pi*x/a)*sin(7*pi*y/b), sin(7*pi*x/a)*sin(1*pi*y/b)];
% dxxphi = diff(phi,x,2); dxyphi = diff(diff(phi,x),y); dyyphi = diff(phi,y,2);
% 
% phi_nm = [  sin(1*pi*x/a)*sin(1*pi*y/b), sin(1*pi*x/a)*sin(3*pi*y/b), ...
%             sin(3*pi*x/a)*sin(1*pi*y/b), sin(3*pi*x/a)*sin(3*pi*y/b), ...
%             sin(1*pi*x/a)*sin(5*pi*y/b), sin(5*pi*x/a)*sin(1*pi*y/b), ...
%             sin(3*pi*x/a)*sin(5*pi*y/b), sin(5*pi*x/a)*sin(3*pi*y/b), ...
%             sin(5*pi*x/a)*sin(5*pi*y/b), sin(1*pi*x/a)*sin(7*pi*y/b), ...
%             sin(7*pi*x/a)*sin(1*pi*y/b), sin(3*pi*x/a)*sin(7*pi*y/b), ...
%             sin(7*pi*x/a)*sin(3*pi*y/b), sin(5*pi*x/a)*sin(7*pi*y/b), ...
%             sin(7*pi*x/a)*sin(5*pi*y/b), sin(7*pi*x/a)*sin(7*pi*y/b)];
% 
% dxxphi_nm = diff(phi_nm,x,2); dxyphi_nm = diff(diff(phi_nm,x),y); dyyphi_nm = diff(phi_nm,y,2);
% 
% % Precompute Global Matrices
% fprintf('Computing Global Matrices...\n');
% Mgm = zeros(Nm); Kgk_undamaged = zeros(Nm); Fgf = zeros(Nm,1);
% for i=1:Nm
%     for j=1:Nm
%         Mgm(i,j) = double(ro*h*int(int(phi(i)*phi(j),x,0,a),y,0,b));
%         Kgk_undamaged(i,j) = double(D*int(int(dxxphi(i)*dxxphi(j)+2*dxyphi(i)*dxyphi(j)+dyyphi(i)*dyyphi(j),x,0,a),y,0,b));
%     end
%     Fgf(i,1) = double(int(int(p*phi(i),x,0,a),y,0,b));
% end
% Cgc = 0.0003*Mgm + 0.0003*Kgk_undamaged;
% 
% % Precompute Damage Influence Matrices (Kd)
% Kd = zeros(Ndam, Nm, Nm);
% for i=1:Ndam
%     for j=1:Nm
%         for k=1:Nm
%             Kd(i,j,k) = double(D*int(int(dxxphi(j)*dxxphi(k)+2*dxyphi(j)*dxyphi(k)+dyyphi(j)*dyyphi(k),x,damdox(i,1),damdox(i,2)),y,damdoy(i,1),damdoy(i,2)));
%         end
%     end
% end
% 
% % M ,K , C for nm time integration
% Mgm_nm = zeros(Nm_nm); Kgk_undamaged_nm = zeros(Nm_nm); Fgf_nm = zeros(Nm_nm,1);
% for i=1:Nm_nm
%     for j=1:Nm_nm
%         Mgm_nm(i,j) = double(ro*h*int(int(phi_nm(i)*phi_nm(j),x,0,a),y,0,b));
%         Kgk_undamaged_nm(i,j) = double(D*int(int(dxxphi_nm(i)*dxxphi_nm(j)+2*dxyphi_nm(i)*dxyphi_nm(j)+dyyphi_nm(i)*dyyphi_nm(j),x,0,a),y,0,b));
%     end
%     Fgf_nm(i,1) = double(int(int(p*phi_nm(i),x,0,a),y,0,b));
% end
% Cgc_nm = 0.0003*Mgm_nm + 0.0003*Kgk_undamaged_nm;
% 
% % Precompute Damage Influence Matrices (Kd)
% Kd_nm = zeros(Ndam, Nm_nm, Nm_nm);
% for i=1:Ndam
%     for j=1:Nm_nm
%         for k=1:Nm_nm
%             Kd_nm(i,j,k) = double(D*int(int(dxxphi_nm(j)*dxxphi_nm(k)+2*dxyphi_nm(j)*dxyphi_nm(k)+dyyphi_nm(j)*dyyphi_nm(k),x,damdox(i,1),damdox(i,2)),y,damdoy(i,1),damdoy(i,2)));
%         end
%     end
% end
% 
% % Sensor Observation Matrix (Phi)
% Phi_mat = zeros(nObserv, Nm);
% for i=1:nObserv
%     x1=sensox(i); y1=sensoy(i);
%     Phi_mat(i,:) = [sin(pi*x1/a)*sin(pi*y1/b), sin(1*pi*x1/a)*sin(3*pi*y1/b), ...
%                     sin(3*pi*x1/a)*sin(1*pi*y1/b), sin(3*pi*x1/a)*sin(3*pi*y1/b), ...
%                     sin(1*pi*x1/a)*sin(5*pi*y1/b), sin(5*pi*x1/a)*sin(1*pi*y1/b), ...
%                     sin(1*pi*x1/a)*sin(7*pi*y1/b), sin(7*pi*x1/a)*sin(1*pi*y1/b)];
% 
% end
% 
% Phi_mat_nm = zeros(nObserv, Nm_nm);
% for i=1:nObserv
%     x1=sensox(i); y1=sensoy(i);
%     Phi_mat_nm(i,:) = [ sin(1*pi*x1/a)*sin(1*pi*y1/b), sin(1*pi*x1/a)*sin(3*pi*y1/b), ...
%                         sin(3*pi*x1/a)*sin(1*pi*y1/b), sin(3*pi*x1/a)*sin(3*pi*y1/b), ...
%                         sin(1*pi*x1/a)*sin(5*pi*y1/b), sin(5*pi*x1/a)*sin(1*pi*y1/b), ...
%                         sin(3*pi*x1/a)*sin(5*pi*y1/b), sin(5*pi*x1/a)*sin(3*pi*y1/b), ...
%                         sin(5*pi*x1/a)*sin(5*pi*y1/b), sin(1*pi*x1/a)*sin(7*pi*y1/b), ...
%                         sin(7*pi*x1/a)*sin(1*pi*y1/b), sin(3*pi*x1/a)*sin(7*pi*y1/b),... 
%                         sin(7*pi*x1/a)*sin(3*pi*y1/b), sin(5*pi*x1/a)*sin(7*pi*y1/b), ...
%                         sin(7*pi*x1/a)*sin(5*pi*y1/b), sin(7*pi*x1/a)*sin(7*pi*y1/b)];
% end
% 
% 
% %% FINDING PARTICIPATION 
% % [V, D] = eig(Kgk_undamaged_nm, Mgm_nm); 
% % [omega2, idx] = sort(diag(D)); % Sort by natural frequency (ascending)
% % V = V(:, idx);   % Reorder mode shapes to match frequencies
% % 
% % omega_sorted = sqrt(omega2);
% % freq_sorted = omega_sorted/(2*pi);
% % % mass norm
% % for i = 1:size(V, 2)
% %     normalization_factor = sqrt(V(:,i)' * Mgm_nm * V(:,i));
% %     V(:,i) = V(:,i) / normalization_factor;
% % end
% % % finding perticipation
% % iota = ones(size(Mgm_nm, 1), 1); % Influence vector
% % Gamma = zeros(size(V, 2), 1);
% % EffMass = zeros(size(V, 2), 1);
% % 
% % for i = 1:size(V, 2)
% %     Gamma(i) = V(:,i)' * Mgm_nm * iota; % Denominator is 1 due to normalization
% %     EffMass(i) = Gamma(i)^2;       % Effective Modal Mass
% % end
% % 
% % % Convert to Percentage
% % ParticipationPercent = (EffMass / sum(EffMass)) * 100;
% 
% %% SAVIG THE SYM EXP
% save('PlateModelData_10NM.mat', 'Mgm', 'Kgk_undamaged', 'Kd', 'Phi_mat','Mgm_nm', 'Kgk_undamaged_nm', 'Kd_nm', 'Phi_mat_nm', 'Fgf', 'Fgf_nm');
% 
% %% LOADING THE DATA 
% load('PlateModelData_10NM.mat');
% %% 2. GENERATE TRUTH DATA (Damaged Simulation)
% % damf_true = zeros(Ndam, N); X_DAM =linspace(0, 0.8, N-199);
% % damf_true(2, 1200:N) = 0.6 * sqrt( (0:N-1200)/(N-1200) );
% % damf_true(2, 1200:N) = linspace(0, 0.6, N-1199);
% % damf_true(5, 30:N) = linspace(0, 0.2, N-29);
% 
% damf_true=zeros(Ndam,N);
% % 6 dam zones 
% % % 
% %   for i=200:N
% %       damf_true(1,i)=(i-200)^.7*0.3/(N-200)^.7;
% %   end
% % % % 
% %   for i=1200:N
% %       damf_true(2,i)=(i-1200)*0.6/(N-1200);
% %   end
% % 
% % %   for i=600:N
% % %      damf_true(3,i)=(i-600)*0.12/(N-600);
% % %  end
% % % 
% % % for i=10:N
% % %     damf_true(4,i)=((i-10)^0.8)*0.3/(N-10);
% % % end
% % 
% % for i=30:N
% %     damf_true(5,i)=(i-30)*0.2/(N-30);
% % end
% % 
% % 
% % for i=300:N
% %     damf_true(6,i)=(i-300)*0.2/(N-300);
% % end
% 
% % for 12 zones 
% for i = 200:N,  damf_true(1,i) = (i-200)*0.3/(N-200); end
% for i = 1200:N, damf_true(2,i) = (i-1200) * 0.1 / (N-1200); end
% for i = 600:N,  damf_true(6,i) = (i-600) * 0.3 / (N-600); end
% for i = 30:N,   damf_true(7,i) = (i-30)*0.2 / (N-30); end
% for i = 800:N,  damf_true(12,i) = (i-800) * 0.1 / (N-800); end
% for i = 1500:N, damf_true(11,i) = (i-1500) * 0.25 / (N-1500); end
% 
% 
% % Run Newmark-Beta to get "Measured" Sensors
% 
% gamma = 0.5; 
% beta = 0.25;
% 
% % Precompute Newmark Constants
% a0 = 1/(beta*deltat^2); a1 = gamma/(beta*deltat);
% a2 = 1/(beta*deltat);   a3 = 1/(2*beta) - 1;
% a4 = gamma/beta - 1;    a5 = deltat/2 * (gamma/beta - 2);
% 
% % Initial Conditions
% x_t = zeros(Nm_nm, 1); 
% v_t = zeros(Nm_nm, 1); 
% % Solve for initial acceleration: a = M \ (F0 - C*v0 - K*x0)
% a_t = Mgm_nm \ (Fgf_nm*sin(w_1*0) - Cgc_nm*v_t - Kgk_undamaged_nm*x_t); 
% zst_measured = zeros(2*nObserv, N);
% zacc_measured = zeros(nObserv, N);
% 
% for i = 2:N
%     t = (i-1)*deltat;
% 
%     % 1. Update current Damage/Stiffness
%     K_curr = Kgk_undamaged_nm;
%     for d = 1:Ndam
%         K_curr = K_curr - damf_true(d,i) * squeeze(Kd_nm(d,:,:)); 
%     end
% 
%     % 2. Effective Stiffness Matrix
%     K_hat = K_curr + a1*Cgc_nm + a0*Mgm_nm;
% 
%     % 3. Effective Force Vector (Predictor)
%     F_ext = Fgf_nm * sin(w_1*t);
%     F_hat = F_ext + Mgm_nm*(a0*x_t + a2*v_t + a3*a_t) + Cgc_nm*(a1*x_t + a4*v_t + a5*a_t);
% 
%     % 4. Solve for New Displacement
%     x_new = K_hat \ F_hat;
% 
%     % 5. Derived Velocity and Acceleration (Corrector)
%     a_new = a0*(x_new - x_t) - a2*v_t - a3*a_t;
%     v_new = v_t + (1-gamma)*deltat*a_t + gamma*deltat*a_new;
% 
%     % 6. Store Observation and update states
%     zst_measured(:,i) = [Phi_mat_nm *( x_new + x_new*0.02*randn); Phi_mat_nm *( v_new + v_new*0.02*randn)];
%     zacc_measured(:,i) = Phi_mat_nm *( a_new + a_new*0.03*randn);
% 
%     x_t = x_new;
%     v_t = v_new;
%     a_t = a_new;
% end
% 
% 
% 
% 
% 
% % storing variables for the DEKF
% 
% 
% 
% 
% 
% 
% 
% %% 3. BAYESIAN OPTIMIZATION
% 
% optimVars_ekf = [%-23.784         -27.9         -20.692        -4.4305      5.9995(0.055485)
% %-16.647        -18.11         -11.75         -7.2362      0.44239(0.044954)
% % -25.333        -26.894        -19.881        -6.139       4.1736(0.052563)
% %-24.654        -27.012        -20.68         -6.6153      4.19833(0.051902)
%     optimizableVariable('p_state_exp', [-26, -22], 'Type', 'real')
%     optimizableVariable('p_param_exp', [-28, -26], 'Type', 'real')
%     optimizableVariable('q_state_exp', [-21, -19], 'Type', 'real')
%     optimizableVariable('q_param_exp', [-8, -4], 'Type', 'real')
%     %optimizableVariable('p_exp', [-10, -6], 'Type', 'real')
%     %optimizableVariable('q_exp', [-20, -12], 'Type', 'real')
%     optimizableVariable('r_exp', [4, 6], 'Type', 'real')
% ];
% 
% ObjFcn_ekf = @(p) run_ekf_internal(p, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                nObserv, Nm, Ndam, deltat, N, w_1, zacc_measured, damf_true);
% 
% 
% 
% optimVars_dkf2 = [%   -12.872 |       -29.96 |       6.0019 |       2.6404 |      -3.1929 |       13.863 |      0.68604 
%     % FOR 12 nb modes 8 sens 
%  % -12.393         -29.406          6.7551          2.9563         -3.8204  13.899       0.20024(0.041737)***
%  % -12.933          -29.58          6.1487          2.4279         -3.8502          13.139       0.8012(0.047118)
% % -13.74          -27.21         -2.5021         -13.005         -3.7231         -1.3997       0.57543(0.051169)
% 
%     optimizableVariable('p_state_dkf2', [-25, -22], 'Type', 'real')
%     optimizableVariable('p_param_dkf2', [-35, -28], 'Type', 'real')
%     optimizableVariable('q_state_dkf2', [2, 8], 'Type', 'real')
%     optimizableVariable('q_param_dkf2', [1, 6], 'Type', 'real')
%     optimizableVariable('r_state_dkf2', [-6, -1], 'Type', 'real')
%     optimizableVariable('r_param_dkf2', [10, 18], 'Type', 'real')
%     optimizableVariable('rval', [0.2, 0.9], 'Type', 'real')
% 
%     ];
%  ObjFcn_dkf2 = @(p) run_dkf2_internal(p, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                         nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);
% 
% 
% optimVars_dkf1 = [% -24.609         -27.247          14.064         -12.041          7.1334         0.99475   (0.061978
%     % FOR 12 nb mode 8 sen 
% % %-25.534         -26.823          14.98          -12.619          5.1501         -1.6602(0.041977)
% % -24.857          -26.32          14.084          -12.05          7.2985         -1.0863(0.04207) 
% % % -25.805         -27.461          15.064         -12.026          8.2364     -0.90961  (0.03678)***
% %
% 
%     optimizableVariable('p_state_dkf1', [-26, -24], 'Type', 'real')
%     optimizableVariable('p_param_dkf1', [-28, -26], 'Type', 'real')
%     optimizableVariable('q_state_dkf1', [14, 15], 'Type', 'real')
%     optimizableVariable('q_param_dkf1', [-14, -12], 'Type', 'real')
%     optimizableVariable('r_state_dkf1', [5, 8], 'Type', 'real')
%     optimizableVariable('r_param_dkf1', [-2, 1], 'Type', 'real')
% 
%     ];
%  ObjFcn_dkf1 = @(p) run_dkf1_internal(p, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                         nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);
% 
% 
% fprintf('\nStarting Bayesian Optimization ekf ...\n');
% results_ekf = bayesopt(ObjFcn_ekf, optimVars_ekf, 'MaxObjectiveEvaluations', 20, 'PlotFcn', {});
% 
% fprintf('\nStarting Bayesian Optimization dkf1 ...\n');
% results_dkf1 = bayesopt(ObjFcn_dkf1, optimVars_dkf1, 'MaxObjectiveEvaluations', 20, 'PlotFcn', {});
% 
% fprintf('\nStarting Bayesian Optimization dkf2 ...\n');
% results_dkf2 = bayesopt(ObjFcn_dkf2, optimVars_dkf2, 'MaxObjectiveEvaluations', 20, 'PlotFcn', {});
% 
% %% 4. FINAL SIMULATION WITH BEST PARAMETERS
% bestP = bestPoint(results_ekf);
% [final_rmse, fd_est ] = run_ekf_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                         nObserv, Nm, Ndam, deltat, N, w_1, zacc_measured, damf_true);
% 
% bestP = bestPoint(results_dkf1);
% % for scheme 2
% [final_rmse_dkf1, fd_dkf1_optimized ] = run_dkf1_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                         nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);
% 
% bestP = bestPoint(results_dkf2);
% % for scheme 2
% [final_rmse_dkf2, fd_dkf2_optimized ] = run_dkf2_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                         nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);
% 
% %%
% figure(1);
% for i=1:Ndam
%     subplot(3,4,i);
%     plot(damf_true(i,:), 'k', 'LineWidth', 1.5); hold on;
%     %plot(fd_dkf2_optimized(i,:), 'b--', 'LineWidth', 1);
%     plot(fd_est(i,:), 'g--', 'LineWidth', 1);
%     title(['Zone ', num2str(i)]); grid on;
%     if i==1; legend('True','EKF Tuned'); end
% end
% %%
% figure(17);
% clf; 
% set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1, 1, 12, 9]); 
% 
% marker_interval = 450; 
% smooth_span = 100;
% smooth_span_ekf = 100;
% for i = 1:Ndam
%     % --- 1. Calculate RMSE ---
%     rmse1 = sqrt(mean((damf_true(i,:) - fd_dkf1_optimized(i,:)).^2));
%     rmse2 = sqrt(mean((damf_true(i,:) - fd_dkf2_optimized(i,:)).^2));
%     rmse3 = sqrt(mean((damf_true(i,:) - fd_est(i,:)).^2));
% 
%     subplot(3,4,i);
%     hold on;
% 
%     % --- 2. Plotting with "Haloed" Markers ---
%     h_true = plot(damf_true(i,:), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 2.5);
% 
%     % DKF1
%     smooth_dk1 = movmean(fd_dkf1_optimized(i,:), smooth_span);
%     h_dk1 = plot(smooth_dk1, 'k-', 'LineWidth', 1.1, ...
%         'Marker', 'o', 'MarkerIndices', 1 : marker_interval : length(smooth_dk1), ...
%         'MarkerSize', 4.5, 'MarkerFaceColor', 'w');
% 
%     % DKF2
%     smooth_dk2 = movmean(fd_dkf2_optimized(i,:), smooth_span);
%     h_dk2 = plot(smooth_dk2, 'k-', 'LineWidth', 1.1, ...
%         'Marker', 's', 'MarkerIndices', floor(marker_interval/3) : marker_interval : length(smooth_dk2), ...
%         'MarkerSize', 4.5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
% 
%     % EKF
%     smooth_ekf = movmean(fd_est(i,:), smooth_span_ekf);
%     h_ekf = plot(smooth_ekf, 'k:.', 'LineWidth', 1.1, ...
%         'Marker', '+', 'MarkerIndices', floor(2*marker_interval/3) : marker_interval : length(smooth_ekf), ...
%         'MarkerSize', 6);
% 
%     % --- 3. Refined RMSE Text Placement ---
%     % Placed at top-left (0.05, 0.9) of each subplot
%     text_str = {sprintf('DK1: %.4f', rmse1), ...
%                 sprintf('DK2: %.4f', rmse2), ...
%                 sprintf('EKF: %.4f', rmse3)};
%     text(0.05, 0.90, text_str, 'Units', 'normalized', ...
%         'VerticalAlignment', 'top', 'FontSize', 7.5, ...
%         'FontName', 'Helvetica', 'BackgroundColor', [1 1 1 0.7]);
% 
%     % --- 4. Subplot Formatting ---
%     title(['Zone ', num2str(i)], 'FontSize', 10, 'FontWeight', 'bold');
%     grid on; box on;
%     set(gca, 'GridAlpha', 0.1, 'FontSize', 9);
%     if mod(i-1, 4) == 0, ylabel('Damage Factor'); end
%     if i > 8, xlabel('Time Step'); end
% 
%     if i == 1, plot_handles = [h_true, h_dk1, h_dk2, h_ekf]; end
% end
% 
% % --- 5. Corrected Legend Placement ---
% % [left, bottom, width, height] in figure-normalized units
% % Values around 0.5 center the legend horizontally
% leg_labels = {'True', 'DKF1[RMSE:0.04207]', 'DKF2[RMSE:0.041737]', ...
%               'EKF[RMSE: 0.052563]  \bf[Sensor:8] (12MODES NM)\rm'};
% L = legend(plot_handles, leg_labels, 'Orientation', 'horizontal', ...
%     'Interpreter', 'tex', 'FontSize', 11);
% set(L, 'Position', [0.15, 0.015, 0.7, 0.03], 'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
% 
% %% --- HELPER: EKF FUNCTION ---
% function [score, fd_history] = run_ekf_internal(p, M, K0, C, Kd_all, Phi, Fgf, nObs, Nm, Ndam, dt, N, w1, z, truth)
%     % Unpack
%     P_state = 10^p.p_state_exp;
%     P_param = 10^p.p_param_exp;
%     Q_state = 10^p.q_state_exp;
%     Q_param = 10^p.q_param_exp;
%     %P_val = 10^(p.p_exp);
%     %Q_val = 10^(p.q_exp);
%     R_val = 10^p.r_exp;
%     state_ekf_store = zeros(Nm,N);
% 
%     % Init
%     xx = zeros(2*Nm + Ndam, 1);
%     %P = blkdiag(1e-6*eye(2*Nm), 1e-6*eye(Ndam));
%     P = 1*blkdiag(P_state*eye(2*Nm), P_param*eye(Ndam));
%     Q = 1*blkdiag(Q_state*eye(2*Nm), Q_param*eye(Ndam));
%     R = R_val * eye(nObs);
%     fd_history = zeros(Ndam, N);
%     Obsv_ekf= zeros(2*Nm+Ndam,1);
% 
% 
% 
% 
% 
%     for i=2:N
%         u = xx(1:Nm); v = xx(Nm+1:2*Nm); fd = xx(2*Nm+1:end);
% 
%         % Jacobian Prep
%         K_eff = K0; S = zeros(Nm, Ndam);
%         for d=1:Ndam
%             Ki = squeeze(Kd_all(d,:,:));
%             K_eff = K_eff - fd(d)*Ki;
%             S(:,d) = Ki * u;
%         end
% 
%         % Continuous Jacobian A
%         A = [zeros(Nm), eye(Nm), zeros(Nm,Ndam);
%              -M\K_eff, -M\C, M\S;
%              zeros(Ndam, 2*Nm+Ndam)];
% 
%         % 2nd Order Taylor Discretization: F = I + Adt + 0.5(Adt)^2
%         At = A*dt;
%         F_jac = eye(2*Nm+Ndam) + At + 0.5*(At*At);
% 
%         % Prediction
%         acc = M \ (Fgf*sin(w1*(i-1)*dt) - K_eff*u - C*v);
%         xx(1:Nm) = u + v*dt;
%         xx(Nm+1:2*Nm) = v + acc*dt;
%         P = F_jac * P * F_jac' + Q;
% 
%         % Update
%         H = [-Phi*(M\K_eff), -Phi*(M\C), Phi*(M\S)];
%         innov = z(:,i) - Phi*acc;
%         K_gain = (P*H') / (H*P*H' + R);
%         xx = xx + K_gain * innov;
%         P = (eye(length(xx)) - K_gain*H) * P;
% 
%         state_ekf_store(:,i)=xx(1:Nm);
%         ob_ekf = obsv(F_jac,H);
% 
% 
%          for k=1: 2*Nm+Ndam
% 
%             Obsv_ekf(k) = Obsv_ekf(k) + norm(ob_ekf(:, k))^2;
% 
%         end
% 
% 
% 
%         fd_history(:,i) = xx(2*Nm+1:end);
%     end
%     figure(10);
%     bar(Obsv_ekf(2*Nm+1:2*Nm+Ndam));%/max(Obsv_fd_s1));
%     set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6','Zone 7', 'Zone 8', 'Zone 9', 'Zone 10', 'Zone 11', 'Zone 12'});
%     title('Observability(Fd) over Time for EKF');
% 
%     score = sqrt(mean((truth(:) - fd_history(:)).^2));
% end
% 
% 
% 
% 
% %% DEKF1 loop 
% function[final_rmse_dkf1, fd_dkf1_opt ] = run_dkf1_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                         nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf)
% 
% X=zeros(nObserv,Nm);
% 
% 
% H=[Phi_mat  X ; X  Phi_mat];
% 
% 
% 
% fac=zeros(Ndam,1);
% facc=fac;
% fd_ex = zeros(Ndam,1);
% fd_est_ex = zeros(Ndam,N);
% fd_dkf2 = zeros(Ndam,1);
% fd_est_dkf2 = zeros(Ndam,N);
% 
% 
% xxp=zeros(2*Nm,1);
% xxp2=zeros(2*Nm,1);
% xxp_exfd=zeros(2*Nm+Ndam,1);  % for the extended kalman filter 
% vvp=zeros(Nm,1);
% aap=zeros(Nm,1);
% accp=zeros(Nm,1);
%    acp=zeros(Nm,1);
%    acp_2 = zeros(Nm,1);
%    xp=zeros(Nm,1);
%    vp=zeros(Nm,1);
%    xrip=zeros(Nm,1);
%    vrip=zeros(Nm,1);
%    accpri=zeros(Nm,1);
% xx_true_store = zeros(Nm,N);   
% xx_store_dkf1 = zeros(Nm,N);
% xx_store_dkf2 = zeros(Nm,N);
% xx_store_ekf = zeros(Nm,N);
% %% Creating some obs matrix for the different schemes and the states ( sc1 & sc2 & extened) 
% Obsv_state_s1 = zeros(2*Nm,1);
% Obsv_fd_s2 = zeros(Ndam,1);
% Obsv_fd_s1 = zeros(Ndam,1);
% Obsv_state_s2 = zeros(2*Nm,1);
% Obsv_ekf = zeros(2*Nm + Ndam, 1);
% 
% 
% % params for dekf
% Pdkf2 = eye(2*Nm)*10^(bestP.p_state_dkf1);
% %Pdkf2 = eye(2*Nm)*10^(-2);
% Qdkf2 = eye(2*Nm)*10^(bestP.q_state_dkf1);
% Pf2 = blkdiag(0.1288,0.0468,0.0494,0.1424,0.1454, 0.0150,0.0151,0.1247,0.1187,0.0621,0.0638,0.1456)*10^(bestP.p_param_dkf1);
% %Pf2 = 1e-5*blkdiag(0.1077 , 0.0544 ,0.1292, 0.1234 ,  0.0522,  0.1149);
% Qf2 = eye(Ndam)*10^(bestP.q_param_dkf1);
% Rdkf2 = eye(2*nObserv)*10^(bestP.r_state_dkf1);
% Rf2 = eye(nObserv)*10^(bestP.r_param_dkf1);
% %r = bestP.rval;
% 
% %  P=[ones(Nm,1)*10^(-2);ones(Nm,1)*10^(-2)];   % 10^-10  for in 4 sensors
% %  P=diag(P);
% %  Pdkf2=P;
% % 
% %  Q=eye(2*Nm)*10^-16;  %  % scheme 1 -> eye(dim)/10^-10 for 6 sensor (unchange modes)
% %  Qdkf2=eye(2*Nm)*10^-16; % scheme 2 -> eye(dim)/10^-8 for 6 sensor (unchange modes)
% % 
% % Pf = 1e-5*blkdiag( 0.2135  , 0.1065 , 0.2160, 0.2371, 0.1551 ,  0.2345);% for 6 sensor scheme 1;
% % Pf2 = 1e-5*blkdiag(0.7579 ,  0.7070 ,0.6272, 0.6903 ,  0.7201,  0.6761); % for 6 sensor scheme 2
% % Qf =  1e-6*Pf; %1e-12*blkdiag( 0.4793  , 0.3464 ,0.5248, 0.4994 ,0.3564, 0.5304);
% % Qf2 = 1e-6*Pf2; %    % 1e1*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);
% % 
% % R=1*eye(2*nObserv)*10^(0.34927);  % for state ( disp & vel) of Scheme - 1 for 6 sen ->1*eye(2*nObserv)*10^(-4) works better ( no noise)  
% % Rf=1*eye(nObserv)*10^(0.34927);
% % Rf2=1*eye(nObserv)*10^(0.34927);
% % %Rf = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);
% % 
% % Rdkf2=1*eye(2*nObserv)*10^(0.34927); % for state ( disp & vel) of Scheme - 2
% %Rf2 = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);
% 
% for i =2:N
%     t = (i-1)*deltat;
%     F=Fgf*sin(w_1*t);    
%     F1=F;
%     Mgminv=inv(Mgm);
%     Af=eye(Ndam);
% 
%     Kgks2 = Kgk_undamaged;
%      for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         Kgks2=Kgks2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
%      end
% 
%       Kdx2=zeros(Nm,Ndam);
% 
%     for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         xi(:,j)=tem*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) );
%         Kdx2(:,j)=xi(:,j);
%     end
% 
%     Hdkf2 = Phi_mat*(inv(Mgm)*(Kdx2)) ;
% 
%     a_dkf2 = Hdkf2*fd_dkf2 + Phi_mat*(inv(Mgm)*(F1- Kgk_undamaged*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm)) - Cgc*(xxp2(Nm+1:2*Nm) + deltat*acp_2) ));
% 
% 
%     Pnf2=Af*Pf2*transpose(Af)+Qf2;
%     K_gain_dkf2=Pnf2*transpose(Hdkf2)*inv((Hdkf2*Pnf2*transpose(Hdkf2))+Rf2);
% 
%     %Af=eye(1);
%     reba=100;
%     fd_dkf2=fd_dkf2+K_gain_dkf2*(zacc_measured(:,i) - a_dkf2);
% 
% 
% 
%       Kgdkf2 = Kgk_undamaged;
%      for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         Kgdkf2=Kgdkf2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
%      end
% 
% 
%         Adkf2=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgdkf2  -inv(Mgm)*Cgc];
%         Bdkf2=[zeros(Nm,Nm); inv(Mgm)];
% 
% 
% 
%     Ax2 = eye(2*Nm) + deltat*Adkf2 + 0.5*deltat^2*(Adkf2*Adkf2);
%     Bx2 = deltat*Bdkf2 ;
% 
%     % % How this is defined
%     % Ax2=deltat*Adkf2+eye(2*Nm);
%     % Bx2=deltat*Bdkf2;
%     % Ax2=zeros(2*Nm,2*Nm);
%     % Bx2_=zeros(2*Nm,2*Nm);
%     % Bx2=zeros(2*Nm,Nm);
%     % for jj=0:100
%     %     Ax2=Ax2+(Adkf2*deltat)^jj/factorial(jj);
%     %     Bx2_=Bx2_+(Adkf2*deltat)^jj/factorial(jj+1);
%     % end
% 
%     %Bx2=Bx2_*Bdkf2*deltat;
%         % For Scheme 2 
%     xx2=Ax2*xxp2+Bx2*F;
%     Pndkf2=Ax2*Pdkf2*transpose(Ax2)+Qdkf2;    
%     Kx2=Pndkf2*transpose(H)*inv((H*Pndkf2*transpose(H))+Rdkf2);
%     xx2=xx2+Kx2*(zst_measured(:,i)-H*xx2);   
%     xx_store_dkf2(:,i) = xx2(1:Nm);
%     fd_est_dkf2(:,i) = fd_dkf2;
% 
%     O_s1_st = obsv(Adkf2,H);
%     O_s1_fd = obsv(Af,Hdkf2);
% 
%     acp_2 = (xx2(Nm+1:2*Nm)-xxp2(Nm+1:2*Nm))/deltat;
%     xxp2 = xx2;
%     Pdkf2=(eye(2*Nm)-Kx2*H)*Pndkf2;
%     Pf2=(eye(Ndam)-K_gain_dkf2*Hdkf2)*Pnf2;
% 
% 
% 
%     % updating the obsv terms
%     for z=1:Ndam
%         Obsv_fd_s1(z) = Obsv_fd_s1(z) + norm(O_s1_fd(:, z))^2;
% 
% 
%     end
% 
%     for k=1: 2*Nm
% 
%         Obsv_state_s1(k) = Obsv_state_s1(k) + norm(O_s1_st(:, k))^2;
% 
%     end
% 
% end
% 
% 
% 
% 
% %% Plotting the observation matrix 
% figure(6)
% bar(Obsv_fd_s1);%/max(Obsv_fd_s1));
% set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6','Zone 7', 'Zone 8', 'Zone 9', 'Zone 10', 'Zone 11', 'Zone 12'});
% title('Observability(Fd) over Time for scheme1');
% 
% 
% 
% figure(8)
% bar(Obsv_state_s1);%/max(Obsv_state_s1));
% set(gca, 'XTickLabel', {'m 1', 'm 2', 'm 3', 'm 4', 'm 5', 'm 6', 'm 7','m 8','m9','m10','m11','m12'});
% title(' Observability Strength (modes) over Time for scheme1 ');
% 
% 
% 
% 
% 
% P1_s1 = 1./Obsv_state_s1;
% P2_s1 = 1./Obsv_fd_s1;
% 
% for i=1:Ndam
%     P2_s1(i);
% end
% 
% %P2_s1 = 1./Obsv_fd_s1;
% 
% %P1_s1 = 1./Obsv_state_s1;
% fd_dkf1_opt = fd_est_dkf2;
% final_rmse_dkf1 = sqrt(mean((damf(:) - fd_est_dkf2(:)).^2));
% end
% 
% 
% 
% %% DEKF2 FUNCTION 
% function[final_rmse_dkf2, fd_dkf2_opt ] = run_dkf2_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                         nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf)
% 
% X=zeros(nObserv,Nm);
% 
% 
% H=[Phi_mat  X ; X  Phi_mat];
% 
% 
% 
% fac=zeros(Ndam,1);
% facc=fac;
% fd_ex = zeros(Ndam,1);
% fd_est_ex = zeros(Ndam,N);
% fd_dkf2 = zeros(Ndam,1);
% fd_est_dkf2 = zeros(Ndam,N);
% 
% 
% xxp=zeros(2*Nm,1);
% xxp2=zeros(2*Nm,1);
% xxp_exfd=zeros(2*Nm+Ndam,1);  % for the extended kalman filter 
% vvp=zeros(Nm,1);
% aap=zeros(Nm,1);
% accp=zeros(Nm,1);
%    acp=zeros(Nm,1);
%    acp_2 = zeros(Nm,1);
%    xp=zeros(Nm,1);
%    vp=zeros(Nm,1);
%    xrip=zeros(Nm,1);
%    vrip=zeros(Nm,1);
%    accpri=zeros(Nm,1);
% xx_true_store = zeros(Nm,N);   
% xx_store_dkf = zeros(Nm,N);
% xx_store_dkf2 = zeros(Nm,N);
% xx_store_ekf = zeros(Nm,N);
% %% Creating some obs matrix for the different schemes and the states ( sc1 & sc2 & extened) 
% Obsv_state_s1 = zeros(2*Nm,1);
% Obsv_fd_s2 = zeros(Ndam,1);
% Obsv_fd_s1 = zeros(Ndam,1);
% Obsv_state_s2 = zeros(2*Nm,1);
% Obsv_ekf = zeros(2*Nm + Ndam, 1);
% 
% 
% % params for dekf
% Pdkf2 = eye(2*Nm)*10^(bestP.p_state_dkf2);
% %Pdkf2 = eye(2*Nm)*10^(-2);
% Qdkf2 = eye(2*Nm)*10^(bestP.q_state_dkf2);
% Pf2 = blkdiag(0.1301,0.04768,0.0502,0.1438,0.1484, 0.0154,0.0155,0.1277,0.1197,0.0631,0.0649,0.1470)*10^(bestP.p_param_dkf2);
% %Pf2 = 1e-5*blkdiag(0.1077 , 0.0544 ,0.1292, 0.1234 ,  0.0522,  0.1149);
% Qf2 = eye(Ndam)*10^(bestP.q_param_dkf2);
% Rdkf2 = eye(2*nObserv)*10^(bestP.r_state_dkf2);
% Rf2 = eye(nObserv)*10^(bestP.r_param_dkf2);
% r = bestP.rval;
% 
% %  P=[ones(Nm,1)*10^(-2);ones(Nm,1)*10^(-2)];   % 10^-10  for in 4 sensors
% %  P=diag(P);
% %  Pdkf2=P;
% % 
% %  Q=eye(2*Nm)*10^-16;  %  % scheme 1 -> eye(dim)/10^-10 for 6 sensor (unchange modes)
% %  Qdkf2=eye(2*Nm)*10^-16; % scheme 2 -> eye(dim)/10^-8 for 6 sensor (unchange modes)
% % 
% % Pf = 1e-5*blkdiag( 0.2135  , 0.1065 , 0.2160, 0.2371, 0.1551 ,  0.2345);% for 6 sensor scheme 1;
% % Pf2 = 1e-5*blkdiag(0.7579 ,  0.7070 ,0.6272, 0.6903 ,  0.7201,  0.6761); % for 6 sensor scheme 2
% % Qf =  1e-6*Pf; %1e-12*blkdiag( 0.4793  , 0.3464 ,0.5248, 0.4994 ,0.3564, 0.5304);
% % Qf2 = 1e-6*Pf2; %    % 1e1*blkdiag(0.0999,0.2832,0.7140 ,0.1923,0.2298,0.2095);
% % 
% % R=1*eye(2*nObserv)*10^(0.34927);  % for state ( disp & vel) of Scheme - 1 for 6 sen ->1*eye(2*nObserv)*10^(-4) works better ( no noise)  
% % Rf=1*eye(nObserv)*10^(0.34927);
% % Rf2=1*eye(nObserv)*10^(0.34927);
% % %Rf = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);
% % 
% % Rdkf2=1*eye(2*nObserv)*10^(0.34927); % for state ( disp & vel) of Scheme - 2
% %Rf2 = 1e-2*blkdiag(0.2524, 0.4019,0.8427,0.4211);
% 
% for i =2:N
%     t = (i-1)*deltat;
%     F=Fgf*sin(w_1*t);    
%     F1=F;
%     Mgminv=inv(Mgm);
%     Af=eye(Ndam);
% 
%     Kgks2 = Kgk_undamaged;
%      for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         Kgks2=Kgks2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
%      end
% 
%       Kdx2=zeros(Nm,Ndam);
% 
%     for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         xi(:,j)=tem*(xxp2(1:Nm) + deltat*xxp2(Nm+1:2*Nm) + deltat^2*r*(1-r)*acp_2 );
%         Kdx2(:,j)=xi(:,j);
%     end
% 
%     Hdkf2 = Phi_mat*(inv(Mgm + deltat*(1-r)*Cgc + deltat^2*(1-r)^2*(Kgks2))*(Kdx2)) ;
% 
%     a_dkf2 = Hdkf2*fd_dkf2 + Phi_mat*(inv(Mgm + deltat*(1-r)*Cgc + deltat^2*(1-r)^2*(Kgks2))*(F1- Kgk_undamaged*(xxp2(1:Nm) + deltat*r*xxp2(Nm+1:2*Nm) + deltat^2*r*(1-r)*acp_2) - Cgc*(xxp2(Nm+1:2*Nm) + deltat*r*acp_2) ));
% 
% 
%     Pnf2=Af*Pf2*transpose(Af)+Qf2;
%     K_gain_dkf2=Pnf2*transpose(Hdkf2)*inv((Hdkf2*Pnf2*transpose(Hdkf2))+Rf2);
% 
%     %Af=eye(1);
%     reba=100;
%     fd_dkf2=fd_dkf2+K_gain_dkf2*(zacc_measured(:,i) - a_dkf2);
% 
% 
% 
%       Kgdkf2 = Kgk_undamaged;
%      for j=1:Ndam
%         for k=1:Nm
%             for m=1:Nm
%                 tem(k,m)=Kd(j,k,m);
%             end
%         end
%         Kgdkf2=Kgdkf2-fd_dkf2(j)*tem;  % here must be fac or facc damf will not be there
%      end
% 
% 
%         Adkf2=[zeros(Nm,Nm)  eye(Nm,Nm);-inv(Mgm)*Kgdkf2  -inv(Mgm)*Cgc];
%         Bdkf2=[zeros(Nm,Nm); inv(Mgm)];
% 
% 
% 
%     Ax2 = eye(2*Nm) + deltat*Adkf2 + 0.5*deltat^2*(Adkf2*Adkf2);
%     Bx2 = deltat*Bdkf2 ;
% 
%     % % How this is defined
%     % Ax2=deltat*Adkf2+eye(2*Nm);
%     % Bx2=deltat*Bdkf2;
%     % Ax2=zeros(2*Nm,2*Nm);
%     % Bx2_=zeros(2*Nm,2*Nm);
%     % Bx2=zeros(2*Nm,Nm);
%     % for jj=0:100
%     %     Ax2=Ax2+(Adkf2*deltat)^jj/factorial(jj);
%     %     Bx2_=Bx2_+(Adkf2*deltat)^jj/factorial(jj+1);
%     % end
% 
%     %Bx2=Bx2_*Bdkf2*deltat;
%         % For Scheme 2 
%     xx2=Ax2*xxp2+Bx2*F;
%     Pndkf2=Ax2*Pdkf2*transpose(Ax2)+Qdkf2;    
%     Kx2=Pndkf2*transpose(H)*inv((H*Pndkf2*transpose(H))+Rdkf2);
%     xx2=xx2+Kx2*(zst_measured(:,i)-H*xx2);   
%     xx_store_dkf2(:,i) = xx2(1:Nm);
%     fd_est_dkf2(:,i) = fd_dkf2;
% 
% 
%     O_s2_st = obsv(Adkf2,H);
%     O_s2_fd = obsv(Af,Hdkf2);
% 
%     acp_2 = (xx2(Nm+1:2*Nm)-xxp2(Nm+1:2*Nm))/deltat;
%     xxp2 = xx2;
%     Pdkf2=(eye(2*Nm)-Kx2*H)*Pndkf2;
%     Pf2=(eye(Ndam)-K_gain_dkf2*Hdkf2)*Pnf2;
%     % updating the obsv terms
%     for z=1:Ndam
%         %Obsv_fd_s1(z) = Obsv_fd_s1(z) + norm(O_s1_fd(:, z))^2;
% 
%         Obsv_fd_s2(z) = Obsv_fd_s2(z) + norm(O_s2_fd(:, z))^2;
%     end
% 
%     for k=1: 2*Nm
% 
%         %Obsv_state_s1(k) = Obsv_state_s1(k) + norm(O_s1_st(:, k))^2;
%         Obsv_state_s2(k) = Obsv_state_s2(k) + norm(O_s2_st(:, k))^2;
%     end
% 
% end
% 
% 
% 
% 
% %% Plotting the observation matrix 
% % figure(6)
% % bar(Obsv_fd_s1);%/max(Obsv_fd_s1));
% % set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6'});
% % title('Observability(Fd) over Time for scheme1');
% 
% figure(7)
% bar(Obsv_fd_s2);%/max(Obsv_fd_s2));
% set(gca, 'XTickLabel', {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5', 'Zone 6','Zone 7', 'Zone 8', 'Zone 9', 'Zone 10', 'Zone 11', 'Zone 12'});
% title('Observability(Fd) over Time for scheme2 ');
% 
% % figure(8)
% % bar(Obsv_state_s1);%/max(Obsv_state_s1));
% % set(gca, 'XTickLabel', {'m 1', 'm 2', 'm 3', 'm 4', 'm 5', 'm 6', 'm 7','m 8','m9','m10','m11','m12'});
% % title(' Observability Strength (modes) over Time for scheme1 ');
% 
% figure(9)
% bar(Obsv_state_s2);%/max(Obsv_state_s2));
% set(gca, 'XTickLabel', {'m 1', 'm 2', 'm 3', 'm 4', 'm 5', 'm 6', 'm 7','m 8','m9','m10','m11','m12'});
% title(' Observability Strength (modes) over Time for scheme2 ');
% 
% 
% 
% P1_s2 = 1./Obsv_state_s2;
% P2_s2 = 1./Obsv_fd_s2;
% 
% for i=1:Ndam
%     P2_s2(i);
% end
% 
% %P2_s1 = 1./Obsv_fd_s1;
% 
% %P1_s1 = 1./Obsv_state_s1;
% fd_dkf2_opt = fd_est_dkf2;
% final_rmse_dkf2 = sqrt(mean((damf(:) - fd_est_dkf2(:)).^2));
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
clear all; close all; clc;
% wORKING BUT HAVE TO CHECK FOR THE DKF2 IS WORKING BETTER OR NOT 

%% 1. PHYSICAL PARAMETERS & SYSTEM SETUP
p=1000; ro=2700; yo=70*10^9; neu=0.3; h=0.001;
a=0.6; b=0.4; Nm=6; Nm_nm=10; Ndam=12; N=5000; deltat=0.0001; w_1=100*pi;
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


%% FINDING MOST OBS LOCATION FOR THE SENSOR 
a = 0.6; b = 0.4;
num_sensors = 12;

% Your specific 6 modes (Indices m, n)
m_list = [1, 1, 2, 1, 3, 3];
n_list = [1, 2, 3, 3, 1, 3]; 
Nm = length(m_list);

% 2. Create a Candidate Grid (Fine grid to find best spots)
% We avoid the exact edges (0 and 1) where displacement is zero
nx = 50; ny = 50;
[X, Y] = meshgrid(linspace(0.05*a, 0.95*a, nx), linspace(0.05*b, 0.95*b, ny));
cand_x = X(:); cand_y = Y(:);
n_cand = length(cand_x);

% 3. Build the Master Modal Matrix (Full Grid)
Psi_full = zeros(n_cand, Nm);
for j = 1:Nm
    Psi_full(:,j) = sin(m_list(j)*pi*cand_x/a) .* sin(n_list(j)*pi*cand_y/b);
end

% 4. Iterative Effective Independence (EfI) Algorithm
% We start with all points and delete the "least informative" one by one
selected_indices = (1:n_cand)';

while length(selected_indices) > num_sensors
    Phi = Psi_full(selected_indices, :);

    % Compute the Orthogonal Projection Matrix (H-matrix)
    % This identifies the contribution of each sensor to the target modes
    E = Phi * ((Phi' * Phi) \ Phi');

    % Find the point with the minimum diagonal value (lowest contribution)
    [~, min_idx] = min(diag(E));

    % Remove that point from the selected set
    selected_indices(min_idx) = [];
end

% 5. Final Results
opt_x = cand_x(selected_indices);
opt_y = cand_y(selected_indices);
Phi_final = Psi_full(selected_indices, :);

% 6. Display Diagnostics
fprintf('--- OPTIMAL PLACEMENT RESULTS ---\n');
fprintf('Final Rank: %d (MUST BE 6)\n', rank(Phi_final));
fprintf('Final Condition Number: %.4f\n', cond(Phi_final));
fprintf('\nFinal Coordinates (x, y):\n');
disp([opt_x, opt_y]);

% 7. Visualization
figure;
plot(opt_x, opt_y, 'rs', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on; grid on;
xlabel('x (m)'); ylabel('y (m)');
title('12 Optimized Asymmetric Sensors for 6 Modes');
axis([0 a 0 b]);

% %% SPECIFICALLY THIS CONFIGURATION 
% % Asymmetric Sensor Placement Coordinates (a=0.6, b=0.4)
% x_sens = [0.1292, 0.1402, 0.1402, 0.2945, 0.2945, 0.2945, 0.3055, 0.3055, 0.3055, 0.4598, 0.4708, 0.4708];
% y_sens = [0.2037, 0.0788, 0.3212, 0.0861, 0.1963, 0.3139, 0.0861, 0.2037, 0.3139, 0.0788, 0.1963, 0.3212];
% 
% % Mode Indices
% m = [1, 1, 2, 1, 3, 3];
% n = [1, 2, 1, 3, 1, 3];
% 
% Phi = zeros(12,6);
% for i = 1:8
%     for j = 1:6
%         Phi(i,j) = sin(m(j)*pi*x_sens(i)/a) * sin(n(j)*pi*y_sens(i)/b);
%     end
% end
% 
% fprintf('Rank: %d\n', rank(Phi)); 
% fprintf('Condition Number: %.4f\n', cond(Phi));
%%
sensox = [0.1292, 0.1402, 0.1402, 0.2945, 0.2945, 0.2945, 0.3055, 0.3055, 0.3055, 0.4598, 0.4708, 0.4708];
sensoy = [0.2037, 0.0788, 0.3212, 0.0861, 0.1963, 0.3139, 0.0861, 0.2037, 0.3139, 0.0788, 0.1963, 0.3212];


 % % sensox=[a/6  5*a/6  a/6  5*a/6   a/2 a/2];
 % % sensoy=[b/4  b/4  3*b/4  3*b/4  b/2  b/4];

%%
% SYMBOLIC EXPRESSION 
% % Modal Shape and Differential Operators
syms x y
phi = [sin(pi*x/a)*sin(pi*y/b), sin(pi*x/a)*sin(2*pi*y/b),...
       sin(2*pi*x/a)*sin(pi*y/b), sin(pi*x/a)*sin(3*pi*y/b),...
       sin(3*pi*x/a)*sin(pi*y/b), sin(3*pi*x/a)*sin(3*pi*y/b)]; ...
       % sin(1*pi*x/a)*sin(5*pi*y/b), sin(5*pi*x/a)*sin(1*pi*y/b),...
       % sin(3*pi*x/a)*sin(5*pi*y/b), sin(5*pi*x/a)*sin(3*pi*y/b),...
       % sin(5*pi*x/a)*sin(5*pi*y/b)];

dxxphi = diff(phi,x,2); dxyphi = diff(diff(phi,x),y); dyyphi = diff(phi,y,2);

phi_nm = [  sin(1*pi*x/a)*sin(1*pi*y/b), sin(1*pi*x/a)*sin(2*pi*y/b), ...
            sin(2*pi*x/a)*sin(1*pi*y/b), sin(2*pi*x/a)*sin(2*pi*y/b), ...
            sin(1*pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(1*pi*y/b), ...
            sin(2*pi*x/a)*sin(3*pi*y/b), sin(3*pi*x/a)*sin(2*pi*y/b), ...
            sin(3*pi*x/a)*sin(3*pi*y/b), sin(1*pi*x/a)*sin(4*pi*y/b)]; ...
            %sin(4*pi*x/a)*sin(1*pi*y/b), sin(2*pi*x/a)*sin(4*pi*y/b)]; ... 
            % sin(4*pi*x/a)*sin(2*pi*y/b), sin(3*pi*x/a)*sin(4*pi*y/b), ...
            % sin(4*pi*x/a)*sin(3*pi*y/b), sin(4*pi*x/a)*sin(4*pi*y/b), ...
            % sin(1*pi*x/a)*sin(5*pi*y/b), sin(5*pi*x/a)*sin(1*pi*y/b),...
            % sin(2*pi*x/a)*sin(5*pi*y/b), sin(5*pi*x/a)*sin(2*pi*y/b),...
            % sin(3*pi*x/a)*sin(5*pi*y/b), sin(5*pi*x/a)*sin(3*pi*y/b), ...
            % sin(4*pi*x/a)*sin(5*pi*y/b), sin(5*pi*x/a)*sin(4*pi*y/b),...
            % sin(5*pi*x/a)*sin(5*pi*y/b), sin(1*pi*x/a)*sin(6*pi*y/b)];

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


% Recompute the Observation Matrix
Phi_mat = zeros(nObserv, Nm);
for i = 1:nObserv
    x1 = sensox(i); y1 = sensoy(i);
    Phi_mat(i,:) = [sin(pi*x1/a)*sin(pi*y1/b), ...
                    sin(pi*x1/a)*sin(2*pi*y1/b), ...
                    sin(2*pi*x1/a)*sin(pi*y1/b), ...
                    sin(1*pi*x1/a)*sin(3*pi*y1/b),sin(3*pi*x1/a)*sin(1*pi*y1/b),sin(3*pi*x1/a)*sin(3*pi*y1/b)];
end


% for i=1:nObserv
%     x1=sensox(i); y1=sensoy(i);
%     Phi_mat(i,:) = [ ...
%         sin(pi*x1/a)*sin(pi*y1/b),   sin(pi*x1/a)*sin(3*pi*y1/b), ...
%         sin(3*pi*x1/a)*sin(pi*y1/b), sin(3*pi*x1/a)*sin(3*pi*y1/b)]; ...
%         %sin(1*pi*x1/a)*sin(5*pi*y1/b), sin(5*pi*x1/a)*sin(1*pi*y1/b), ...
%         %sin(3*pi*x1/a)*sin(5*pi*y1/b), sin(5*pi*x1/a)*sin(3*pi*y1/b), ...
%         %sin(5*pi*x1/a)*sin(5*pi*y1/b)];
% end

Phi_mat_nm = zeros(nObserv, Nm_nm);
for i=1:nObserv
    x1=sensox(i); y1=sensoy(i);
      % a, b are plate dimensions
    Phi_mat_nm(i,:) = [ ...
        sin(1*pi*x1/a)*sin(1*pi*y1/b), sin(1*pi*x1/a)*sin(2*pi*y1/b), ...
        sin(2*pi*x1/a)*sin(1*pi*y1/b), sin(2*pi*x1/a)*sin(2*pi*y1/b), ...
        sin(1*pi*x1/a)*sin(3*pi*y1/b), sin(3*pi*x1/a)*sin(1*pi*y1/b), ...
        sin(2*pi*x1/a)*sin(3*pi*y1/b), sin(3*pi*x1/a)*sin(2*pi*y1/b), ...
        sin(3*pi*x1/a)*sin(3*pi*y1/b), sin(1*pi*x1/a)*sin(4*pi*y1/b)]; ...
        %sin(4*pi*x1/a)*sin(1*pi*y1/b), sin(2*pi*x1/a)*sin(4*pi*y1/b)]; ... 
        % sin(4*pi*x1/a)*sin(2*pi*y1/b), sin(3*pi*x1/a)*sin(4*pi*y1/b), ...
        % sin(4*pi*x1/a)*sin(3*pi*y1/b), sin(4*pi*x1/a)*sin(4*pi*y1/b), ...
        % sin(1*pi*x1/a)*sin(5*pi*y1/b), sin(5*pi*x1/a)*sin(1*pi*y1/b), ...
        % sin(2*pi*x1/a)*sin(5*pi*y1/b), sin(5*pi*x1/a)*sin(2*pi*y1/b), ...
        % sin(3*pi*x1/a)*sin(5*pi*y1/b), sin(5*pi*x1/a)*sin(3*pi*y1/b), ...
        % sin(4*pi*x1/a)*sin(5*pi*y1/b), sin(5*pi*x1/a)*sin(4*pi*y1/b), ...
        % sin(5*pi*x1/a)*sin(5*pi*y1/b), sin(1*pi*x1/a)*sin(6*pi*y1/b)];
end

% % =========================================================================
% % ANALYTICAL MODAL DEFINITIONS (No Symbolic Math)
% % Define the mode number pairs (m, n) for the 9-mode and 26-mode sets
% % =========================================================================
% modes_9 = [1,1; 1,2; 2,1; 2,2; 1,3; 3,1; 2,3; 3,2; 3,3; 1,4; ...
%             4,1; 2,4; 4,2; 3,4; 4,3; 4,4; 1,5; 5,1; 2,5; 5,2; ...
%             3,5; 5,3; 4,5; 5,4; 5,5; 1,6];
% Nm = size(modes_9, 1);
% 
% modes_26 = [1,1; 1,2; 2,1; 2,2; 1,3; 3,1; 2,3; 3,2; 3,3; 1,4; ...
%             4,1; 2,4; 4,2; 3,4; 4,3; 4,4; 1,5; 5,1; 2,5; 5,2; ...
%             3,5; 5,3; 4,5; 5,4; 5,5; 1,6];
% Nm_nm = size(modes_26, 1);
% 
% fprintf('Computing Analytical Global Matrices...\n');
% 
% % =========================================================================
% % 1. GLOBAL MATRICES FOR 9 MODES (Nm = 9)
% % =========================================================================
% Mgm = zeros(Nm); Kgk_undamaged = zeros(Nm); Fgf = zeros(Nm,1);
% 
% for i = 1:Nm
%     m1 = modes_9(i,1); n1 = modes_9(i,2);
%     alpha = m1*pi/a; beta = n1*pi/b;
% 
%     for j = 1:Nm
%         m2 = modes_9(j,1); n2 = modes_9(j,2);
% 
%         % Over the entire plate, modes are orthogonal. 
%         % Off-diagonals are exactly zero.
%         if (m1 == m2) && (n1 == n2)
%             Mgm(i,j) = ro * h * (a/2) * (b/2);
%             Kgk_undamaged(i,j) = D * (a/2) * (b/2) * (alpha^2 + beta^2)^2;
%         end
%     end
% 
%     % Analytical integral of Force (sin from 0 to L)
%     % is L/(m*pi)*(1-cos(m*pi)). It is 0 for even m, and 2L/(m*pi) for odd.
%     int_x = (a/(m1*pi)) * (1 - cos(m1*pi));
%     int_y = (b/(n1*pi)) * (1 - cos(n1*pi));
%     Fgf(i,1) = p * int_x * int_y;
% end
% 
% Cgc = 0.0003 * Mgm + 0.0003 * Kgk_undamaged;
% 
% % --- Precompute Damage Influence Matrices (Kd) for 9 modes ---
% Kd = zeros(Ndam, Nm, Nm);
% for idx = 1:Ndam
%     x1 = damdox(idx,1); x2 = damdox(idx,2);
%     y1 = damdoy(idx,1); y2 = damdoy(idx,2);
% 
%     for i = 1:Nm
%         m1 = modes_9(i,1); n1 = modes_9(i,2);
%         a1 = m1*pi/a; b1 = n1*pi/b;
% 
%         for j = 1:Nm
%             m2 = modes_9(j,1); n2 = modes_9(j,2);
%             a2 = m2*pi/a; b2 = n2*pi/b;
% 
%             % Use analytical integration functions (defined at bottom)
%             Sx = int_sin_sin(m1, m2, a, x1, x2);
%             Sy = int_sin_sin(n1, n2, b, y1, y2);
%             Cx = int_cos_cos(m1, m2, a, x1, x2);
%             Cy = int_cos_cos(n1, n2, b, y1, y2);
% 
%             term1 = (a1^2 * a2^2 + b1^2 * b2^2) * Sx * Sy;
%             term2 = 2 * (a1 * b1 * a2 * b2) * Cx * Cy;
% 
%             Kd(idx, i, j) = D * (term1 + term2);
%         end
%     end
% end
% 
% % =========================================================================
% % 2. GLOBAL MATRICES FOR 26 MODES (Nm_nm = 26)
% % =========================================================================
% Mgm_nm = zeros(Nm_nm); Kgk_undamaged_nm = zeros(Nm_nm); Fgf_nm = zeros(Nm_nm,1);
% 
% for i = 1:Nm_nm
%     m1 = modes_26(i,1); n1 = modes_26(i,2);
%     alpha = m1*pi/a; beta = n1*pi/b;
% 
%     for j = 1:Nm_nm
%         m2 = modes_26(j,1); n2 = modes_26(j,2);
%         if (m1 == m2) && (n1 == n2)
%             Mgm_nm(i,j) = ro * h * (a/2) * (b/2);
%             Kgk_undamaged_nm(i,j) = D * (a/2) * (b/2) * (alpha^2 + beta^2)^2;
%         end
%     end
%     int_x = (a/(m1*pi)) * (1 - cos(m1*pi));
%     int_y = (b/(n1*pi)) * (1 - cos(n1*pi));
%     Fgf_nm(i,1) = p * int_x * int_y;
% end
% 
% Cgc_nm = 0.0003 * Mgm_nm + 0.0003 * Kgk_undamaged_nm;
% 
% % --- Precompute Damage Influence Matrices (Kd_nm) for 26 modes ---
% Kd_nm = zeros(Ndam, Nm_nm, Nm_nm);
% for idx = 1:Ndam
%     x1 = damdox(idx,1); x2 = damdox(idx,2);
%     y1 = damdoy(idx,1); y2 = damdoy(idx,2);
% 
%     for i = 1:Nm_nm
%         m1 = modes_26(i,1); n1 = modes_26(i,2);
%         a1 = m1*pi/a; b1 = n1*pi/b;
% 
%         for j = 1:Nm_nm
%             m2 = modes_26(j,1); n2 = modes_26(j,2);
%             a2 = m2*pi/a; b2 = n2*pi/b;
% 
%             Sx = int_sin_sin(m1, m2, a, x1, x2);
%             Sy = int_sin_sin(n1, n2, b, y1, y2);
%             Cx = int_cos_cos(m1, m2, a, x1, x2);
%             Cy = int_cos_cos(n1, n2, b, y1, y2);
% 
%             term1 = (a1^2 * a2^2 + b1^2 * b2^2) * Sx * Sy;
%             term2 = 2 * (a1 * b1 * a2 * b2) * Cx * Cy;
% 
%             Kd_nm(idx, i, j) = D * (term1 + term2);
%         end
%     end
% end
% 
% % =========================================================================
% % 3. SENSOR OBSERVATION MATRICES (PHI)
% % =========================================================================
% % Cleanly generated by looping through the mode index arrays
% Phi_mat = zeros(nObserv, Nm);
% for i = 1:nObserv
%     x_s = sensox(i); y_s = sensoy(i);
%     for j = 1:Nm
%         m = modes_9(j,1); n = modes_9(j,2);
%         Phi_mat(i, j) = sin(m*pi*x_s/a) * sin(n*pi*y_s/b);
%     end
% end
% 
% Phi_mat_nm = zeros(nObserv, Nm_nm);
% for i = 1:nObserv
%     x_s = sensox(i); y_s = sensoy(i);
%     for j = 1:Nm_nm
%         m = modes_26(j,1); n = modes_26(j,2);
%         Phi_mat_nm(i, j) = sin(m*pi*x_s/a) * sin(n*pi*y_s/b);
%     end
% end
% 
% 
% % =========================================================================
% % LOCAL FUNCTIONS FOR ANALYTICAL INTEGRATION
% % (Place these at the very end of your .m file)
% % =========================================================================
% 
% % Evaluates the exact integral of sin(m*pi*x/L)*sin(n*pi*x/L) from x1 to x2
% function I = int_sin_sin(m, n, L, x1, x2)
%     if m == n
%         I = 0.5 * ((x2 - x1) - (L/(2*m*pi)) * (sin(2*m*pi*x2/L) - sin(2*m*pi*x1/L)));
%     else
%         I = 0.5 * ( (L/((m-n)*pi)) * (sin((m-n)*pi*x2/L) - sin((m-n)*pi*x1/L)) ...
%                   - (L/((m+n)*pi)) * (sin((m+n)*pi*x2/L) - sin((m+n)*pi*x1/L)) );
%     end
% end
% 
% % Evaluates the exact integral of cos(m*pi*x/L)*cos(n*pi*x/L) from x1 to x2
% function I = int_cos_cos(m, n, L, x1, x2)
%     if m == n
%         I = 0.5 * ((x2 - x1) + (L/(2*m*pi)) * (sin(2*m*pi*x2/L) - sin(2*m*pi*x1/L)));
%     else
%         I = 0.5 * ( (L/((m-n)*pi)) * (sin((m-n)*pi*x2/L) - sin((m-n)*pi*x1/L)) ...
%                   + (L/((m+n)*pi)) * (sin((m+n)*pi*x2/L) - sin((m+n)*pi*x1/L)) );
%     end
% end

%% FINDING PARTICIPATION 
% 1. Define the 26-mode pairs for mapping
m_list = [1,1,2,2,1,3,2,3,3,1,4,2,4,3,4,4,1,5,2,5,3,5,4,5,5,1];
n_list = [1,2,1,2,3,1,3,2,3,4,1,4,2,4,3,4,5,1,5,2,5,3,5,4,5,6];

% 2. Calculate Analytical Participation (Spatial Integral)
Gamma_analytical = zeros(Nm_nm, 1);
for i = 1:Nm_nm
    m = m_list(i); 
    n = n_list(i);
    
    % The physical work done by a uniform load is the integral of the sine waves
    int_x = (a/(m*pi)) * (1 - cos(m*pi));
    int_y = (b/(n*pi)) * (1 - cos(n*pi));
    
    % Gamma is the spatial integral
    Gamma_analytical(i) = int_x * int_y;
end

% 3. Calculate Percentage based on Effective Mass (Gamma^2)
EffMass = Gamma_analytical.^2;
ParticipationPercent = (EffMass / sum(EffMass)) * 100;

% 4. Identify the indices for your EKF (Top 4 Odd-Odd Modes)
[sortedP, sortIdx] = sort(ParticipationPercent, 'descend');

fprintf('\n--- PHYSICAL MODAL PARTICIPATION ---\n');
fprintf('Rank | Mode (m,n) | Participation %% | EKF Priority\n');
for i = 1:6
    idx = sortIdx(i);
    if sortedP(i) > 0.01  % Filter out the 0% (Even) modes
        fprintf('%4d | (%d,%d)      | %10.2f%%    | HIGH\n', ...
            i, m_list(idx), n_list(idx), sortedP(i));
    end
end
%% SAVIG THE SYM EXP
save('PlateModelData_10NM.mat', 'Mgm', 'Kgk_undamaged', 'Kd', 'Phi_mat','Mgm_nm', 'Kgk_undamaged_nm', 'Kd_nm', 'Phi_mat_nm', 'Fgf', 'Fgf_nm');

%% LOADING THE DATA 
load('PlateModelData_10NM.mat');
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
for i = 200:N,  damf_true(1,i) = (i-200)^0.99*0.3/(N-200)^0.99; end
for i = 1200:N, damf_true(2,i) = (i-1200) * 0.1 / (N-1200); end
for i = 600:N,  damf_true(6,i) = (i-600) * 0.3 / (N-600); end
for i = 30:N,   damf_true(7,i) = (i-30)*0.2 / (N-30); end
for i = 800:N,  damf_true(11,i) = (i-800) * 0.1 / (N-800); end
for i = 1500:N, damf_true(10,i) = (i-1500) * 0.25 / (N-1500); end


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
%pSO
% 0.77621        -14.684        -27.009        -11.541       -1  (0.0573)
% 8.1134         -14.863        -24.724        -9.0787      1.5397(0.0571)

% 12nm and 6kf MODES     
%-12.854        -14.092        -28.061        -23.763      -7.2308(0.050954)
%-17.628        -15.991        -23.689        -7.3843      3.3391(0.057785)


%-16.647        -18.11         -11.75         -7.2362      0.44239(0.044954)
% -25.333        -26.894        -19.881        -6.139       4.1736(0.052563)
%-24.654        -27.012        -20.68         -6.6153      4.19833(0.051902)
    optimizableVariable('p_state_exp', [5, 10], 'Type', 'real')
    optimizableVariable('p_param_exp', [-15, -12], 'Type', 'real')
    optimizableVariable('q_state_exp', [-26, -20], 'Type', 'real')
    optimizableVariable('q_param_exp', [-12, -8], 'Type', 'real')
    %optimizableVariable('p_exp', [-10, -6], 'Type', 'real')
    %optimizableVariable('q_exp', [-20, -12], 'Type', 'real')
    optimizableVariable('r_exp', [-1, 4], 'Type', 'real')
];




optimVars_dkf2 = [
    %pso
%-14           -17.829            -2              8            -2.6173            16         0.9(0.0332) 
% -24           -20.859         -1.3989          8.5075            -2            16.402       0.9 (0.0328)
    % BAYOPT
    % -13.664         -28.637          2.7683          7.9198         -0.98527         15.033       0.84917(0.043712)
    %-12.003         -24.921         0.77532          7.8358         -2.0141          14.921       0.86198(0.041496)
 

    optimizableVariable('p_state_dkf2', [-28, -24], 'Type', 'real')
    optimizableVariable('p_param_dkf2', [-22, -18], 'Type', 'real')
    optimizableVariable('q_state_dkf2', [-4, 0], 'Type', 'real')
    optimizableVariable('q_param_dkf2', [6, 10], 'Type', 'real')
    optimizableVariable('r_state_dkf2', [-2, 2], 'Type', 'real')
    optimizableVariable('r_param_dkf2', [16, 20], 'Type', 'real')
    optimizableVariable('rval', [0.4, 0.9], 'Type', 'real')

    ];
 


optimVars_dkf1 = [
    % BAYOPT
%   -28.544         -15.695         -0.37309        -3.8688         -2.4984          2.9227  (0.044715)
% -25.116         -16.068          -1.911            -6            -3.86           2.0535(0.026599)

%% pSO 
%-37.831         -13.449          1.1731          -13.77         0.47721         -6.1478 (0.0283) 
% -36             -15           -1.5629           -10            -2.242            -2 (0.0329)
    % FOR 12 nb mode 8 sen 
% %-25.534         -26.823          14.98          -12.619          5.1501         -1.6602(0.041977)
% -25.116         -16.068          -1.911            -6            -3.86           2.0535(0.026599) 
% % -25.805         -27.461          15.064         -12.026          8.2364     -0.90961  (0.03678)***


    optimizableVariable('p_state_dkf1', [-39, -25], 'Type', 'real')
    optimizableVariable('p_param_dkf1', [-15, -11], 'Type', 'real')
    optimizableVariable('q_state_dkf1', [-1, 2], 'Type', 'real')
    optimizableVariable('q_param_dkf1', [-16, -12], 'Type', 'real')
    optimizableVariable('r_state_dkf1', [-1, 2], 'Type', 'real')
    optimizableVariable('r_param_dkf1', [-2,2], 'Type', 'real')
    
    ];

ObjFcn_ekf = @(p) run_ekf_internal(p, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                               nObserv, Nm, Ndam, deltat, N, w_1, zacc_measured, damf_true);


 ObjFcn_dkf1 = @(p) run_dkf1_internal(p, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);
 ObjFcn_dkf2 = @(p) run_dkf2_internal(p, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);
%% for the bayopt 



% fprintf('\nStarting Bayesian Optimization ekf ...\n');
% results_ekf = bayesopt(ObjFcn_ekf, optimVars_ekf, 'MaxObjectiveEvaluations', 300, 'PlotFcn', {});
% 
% fprintf('\nStarting Bayesian Optimization dkf1 ...\n');
% results_dkf1 = bayesopt(ObjFcn_dkf1, optimVars_dkf1, 'MaxObjectiveEvaluations', 300, 'PlotFcn', {});
% 
% fprintf('\nStarting Bayesian Optimization dkf2 ...\n');
% results_dkf2 = bayesopt(ObjFcn_dkf2, optimVars_dkf2, 'MaxObjectiveEvaluations', 300, 'PlotFcn', {});

%% for PSO 
psoOptions = optimoptions('particleswarm', ...
    'SwarmSize', 10, ...          % 30 particles flying at once
    'MaxIterations', 10, ...      % 30 * 20 = 600 total evaluations, but much faster due to parallelization
    'Display', 'iter', ...
    'UseParallel', false);         % CRITICAL: Evaluates multiple particles simultaneously on your CPU cores

% 1. Dynamically extract lower bounds (lb), upper bounds (ub), and names from your optimVars
% --- EKF ---
n_ekf = length(optimVars_ekf);
lb_ekf = zeros(1, n_ekf); ub_ekf = zeros(1, n_ekf); names_ekf = cell(1, n_ekf);
for i=1:n_ekf
    lb_ekf(i) = optimVars_ekf(i).Range(1); 
    ub_ekf(i) = optimVars_ekf(i).Range(2); 
    names_ekf{i} = optimVars_ekf(i).Name; 
end

% --- DKF1 ---
n_dkf1 = length(optimVars_dkf1);
lb_dkf1 = zeros(1, n_dkf1); ub_dkf1 = zeros(1, n_dkf1); names_dkf1 = cell(1, n_dkf1);
for i=1:n_dkf1
    lb_dkf1(i) = optimVars_dkf1(i).Range(1); 
    ub_dkf1(i) = optimVars_dkf1(i).Range(2); 
    names_dkf1{i} = optimVars_dkf1(i).Name; 
end

% --- DKF2 ---
n_dkf2 = length(optimVars_dkf2);
lb_dkf2 = zeros(1, n_dkf2); ub_dkf2 = zeros(1, n_dkf2); names_dkf2 = cell(1, n_dkf2);
for i=1:n_dkf2
    lb_dkf2(i) = optimVars_dkf2(i).Range(1); 
    ub_dkf2(i) = optimVars_dkf2(i).Range(2); 
    names_dkf2{i} = optimVars_dkf2(i).Name; 
end

% 2. Create Wrapper Functions to convert PSO's numeric vector back into the Table your functions expect
ObjFcn_ekf_PSO  = @(x) ObjFcn_ekf(array2table(x, 'VariableNames', names_ekf));
ObjFcn_dkf1_PSO = @(x) ObjFcn_dkf1(array2table(x, 'VariableNames', names_dkf1));
ObjFcn_dkf2_PSO = @(x) ObjFcn_dkf2(array2table(x, 'VariableNames', names_dkf2));

% 3. Execute Particle Swarm Optimization
fprintf('\nStarting PSO Optimization for EKF ...\n');
[x_ekf_vec, fval_ekf] = particleswarm(ObjFcn_ekf_PSO, n_ekf, lb_ekf, ub_ekf, psoOptions);
% Convert best result back to table format so your downstream code can use results_ekf normally
results_ekf = array2table(x_ekf_vec, 'VariableNames', names_ekf);

fprintf('\nStarting PSO Optimization for DKF1 ...\n');
[x_dkf1_vec, fval_dkf1] = particleswarm(ObjFcn_dkf1_PSO, n_dkf1, lb_dkf1, ub_dkf1, psoOptions);
results_dkf1 = array2table(x_dkf1_vec, 'VariableNames', names_dkf1);

fprintf('\nStarting PSO Optimization for DKF2 ...\n');
[x_dkf2_vec, fval_dkf2] = particleswarm(ObjFcn_dkf2_PSO, n_dkf2, lb_dkf2, ub_dkf2, psoOptions);
results_dkf2 = array2table(x_dkf2_vec, 'VariableNames', names_dkf2);

fprintf('\nOptimization Complete!\n');

%% 4. FINAL SIMULATION WITH BEST PARAMETERS
%% for PSO 
fprintf('\nEvaluating optimized EKF...\n');
[final_rmse, fd_est] = run_ekf_internal(results_ekf, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1, zacc_measured, damf_true);

% Evaluate optimized DKF Scheme 1
fprintf('\nEvaluating optimized DKF1...\n');
[final_rmse_dkf1, fd_dkf1_optimized] = run_dkf1_internal(results_dkf1, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1, zst_measured, zacc_measured, damf_true);

% Evaluate optimized DKF Scheme 2
fprintf('\nEvaluating optimized DKF2...\n');
[final_rmse_dkf2, fd_dkf2_optimized] = run_dkf2_internal(results_dkf2, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
                                        nObserv, Nm, Ndam, deltat, N, w_1, zst_measured, zacc_measured, damf_true);
disp('--- Best EKF Parameters Found ---');
disp(results_ekf);
fprintf('\nFinal RMSE EKF:  %.4f\n', final_rmse);
disp('--- Best DKF1 Parameters Found ---');
disp(results_dkf1);
fprintf('Final RMSE DKF1: %.4f\n', final_rmse_dkf1);
disp('--- Best DKF2 Parameters Found ---');
disp(results_dkf2);
fprintf('Final RMSE DKF2: %.4f\n', final_rmse_dkf2);


%% bayopt operation 
% bestP = bestPoint(results_ekf);
% [final_rmse, fd_est ] = run_ekf_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                         nObserv, Nm, Ndam, deltat, N, w_1, zacc_measured, damf_true);
% 
% bestP = bestPoint(results_dkf1);
% % for scheme 2
% [final_rmse_dkf1, fd_dkf1_optimized ] = run_dkf1_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                         nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);
% 
% bestP = bestPoint(results_dkf2);
% % for scheme 2
% [final_rmse_dkf2, fd_dkf2_optimized ] = run_dkf2_internal(bestP, Mgm, Kgk_undamaged, Cgc, Kd, Phi_mat, Fgf, ...
%                                         nObserv, Nm, Ndam, deltat, N, w_1,zst_measured, zacc_measured, damf_true);

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
leg_labels = {'True', 'DKF1[RMSE:0.0451]', 'DKF2[RMSE:0.0388]', ...
              'EKF[RMSE: 0.0581]  \bf[Sensor:12] (12MODES NM)\rm'};
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
   %rank(Hdkf2)
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
   % rank(Hdkf2)
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
