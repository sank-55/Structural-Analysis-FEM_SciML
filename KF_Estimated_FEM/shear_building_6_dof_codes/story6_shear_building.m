%% =========================================================
%  6-DOF SHEAR BUILDING
%  Natural Frequency Calculation
%
%  Souvik Sarkar
%  Model  : Lumped Parameter Model
%  
% ==========================================================

clear;
clc;
close all;  
  % time interval of the acquisition system
%% Generate the xlxs file ( There toral 5 set of reading )  
filename = 'SB16SEPUDR2.xlsx';  % undamaged 
%filename = 'SB17SEPDMR2.xlsx';   % damaged 
% Import data: Skip the first 8 lines of metadata headers
opts = detectImportOptions(filename, 'NumHeaderLines', 8);
data = readtable(filename, opts);

% Extract columns (Assuming 1st column is Time, next 6 are Accelerations)
raw_time = data{3:end, 1};              % 1st file 1sts 2 data sets are NaN 
raw_accel = data{3:end, 2:7}; % All 6 signals
raw_time = raw_time - raw_time(1);

%% checking any unwanted values   
% Use different names for the variables
[rows, cols] = find(isnan(raw_accel));

% Using 'TimeIndex' and 'ChannelIndex' to avoid conflicts
nan_locations = table(rows, cols, 'VariableNames', {'TimeIndex', 'ChannelIndex'});

disp('Locations of NaN values:');
disp(nan_locations);


%% ---------------------------------------------------------
% GENERATING THE VELOCITY AND THE DISPLACEMEN MEASUREMNET MODEL OUT OF IT  
% ----------------------------------------------------------
dt = 0.00024;
fs = 1/dt; % Example value, replace with your actual Fs
fc_high = 100; % high cutoff 1.2 -1.5 times highest freq % Cut-off frequency in Hz (adjust based on your signal range)
fc_low = 2.5;  % low Cuts off  0.5 or 0.7 times 1st freq Hz
[b, a] = butter(4, [fc_low,fc_high]/(fs/2), 'bandpass');

% Remove trends and filter each column
accel_detrend = detrend(raw_accel); 
accel_filt = filtfilt(b, a, accel_detrend);

% 2. First Integration (Acceleration to Velocity)
velocity = cumtrapz(raw_time, accel_filt);

% 3. Optional: Filter velocity to further reduce drift
velocity_filt = filtfilt(b, a, velocity);

% 4. Second Integration (Velocity to Displacement)
displacement = cumtrapz(raw_time, velocity_filt);

% % Plotting to verify accuracy
% figure('Color', 'w');
% subplot(3,1,1); plot(raw_time, accel_filt, 'k'); ylabel('Accel (m/s^2)'); title('Acceleration (Cleaned)');
% subplot(3,1,2); plot(raw_time, velocity_filt, 'b'); ylabel('Vel (m/s)'); title('Velocity');
% subplot(3,1,3); plot(raw_time, displacement, 'r'); ylabel('Disp (m)'); title('Displacement');
% xlabel('Time (s)');

%% Plot and check 
%% fft 
nfft = 2048;                   % FFT resolution block size
win = hamming(nfft);           % Hamming window to reduce spectral leakage
noverlap = nfft / 2;           % 50% overlap
accel_filt_=accel_filt';
[num_channels, ~] = size(accel_filt_); % 6 floor channels

% Obtain the frequency axis vector
[~, f_axis] = cpsd(accel_filt_(1,:), accel_filt_( 6,:), win, noverlap, nfft, fs);
num_freqs = length(f_axis);

% Pre-allocate 3D CPSD Matrix (6 x 6 x num_freqs)
Pyy = zeros(num_channels, num_channels, num_freqs);

% Build 6x6 CPSD matrix for every channel pair (i, j)
for i = 1:num_channels
    for j = 1:num_channels
        Pyy(i, j, :) = cpsd(accel_filt_(i, :), accel_filt_(j, :), win, noverlap, nfft, fs);
    end
end

% Perform Singular Value Decomposition (SVD) at Each Frequency Line
S1 = zeros(num_freqs, 1);
U1 = zeros(num_channels, num_freqs);

for k = 1:num_freqs
    [U, S, ~] = svd(Pyy(:, :, k));
    S1(k) = S(1, 1);           % Dominant singular value (energy peak)
    U1(:, k) = U(:, 1);         % Dominant singular vector (mode shape candidate)
end

% Plot Singular Value Spectrum to Identify Natural Frequencies
figure('Name', 'FDD Peak Picking Spectrum', 'Color', 'w');
semilogy(f_axis, S1, 'b-', 'LineWidth', 1.5);
grid on; 
xlabel('Frequency (Hz)'); 
ylabel('1st Singular Value');
title('FDD Peak Picking for Natural Frequencies');
xlim([0, fs/2]);
%% plotting undamaged data   
%ffd_data_undamaged = [f_axis, S1];  % undamaged data is stored
figure('Name', 'FDD Peak Picking Spectrum', 'Color', 'w');
semilogy(ffd_data_undamaged(:,1), ffd_data_undamaged(:,2), f_axis, S1);
grid on; 
xlabel('Frequency (Hz)'); 
ylabel('1st Singular Value');
title('FDD Peak Picking for Natural Frequencies Undamaged Case');
xlim([0, fs/2]);
%%  
%
num_modes = 6; % Number of dominant modes to extract (e.g., 3 lowest modes)

df = f_axis(2) - f_axis(1); 

% 2. Convert S1 to log scale so prominence works uniformly across low and high modes
log_S1 = log10(S1);

% 3. Find peaks based on local prominence rather than absolute height
[~, peak_indices] = findpeaks(log_S1, ...
                              'MinPeakDistance', floor(1.5 / df), ... % At least 1.5 Hz apart
                              'MinPeakProminence', 0.4);               % Stands out by at least 0.4 decades

% 4. Keep only up to num_modes and ensure ascending frequency order
if length(peak_indices) > num_modes
    % Sort identified peaks by prominence to pick top 'num_modes', then re-sort by frequency
    [~, prom_order] = sort(log_S1(peak_indices), 'descend');
    peak_indices = peak_indices(prom_order(1:num_modes));
    peak_indices = sort(peak_indices); % Restore ascending order (Mode 1, Mode 2, ...)
end
natural_freqs = f_axis(peak_indices);
fprintf('\n--- Identified Natural Frequencies ---\n');
for m = 1:length(natural_freqs)
    fprintf('Mode %d: %.2f Hz\n', m, natural_freqs(m));
end

% =========================================================================
% EXTRACT AND NORMALIZE EXP. MODE SHAPES
% =========================================================================
% Pre-allocate mode shape matrix (6 floors x num_modes)
phi_exp = zeros(num_channels, length(natural_freqs));

for m = 1:length(natural_freqs)
    idx = peak_indices(m);
    
    % Extract complex singular vector at peak frequency
    u_complex = U1(:, idx);
    
    % FDD singular vectors can have a complex phase angle due to noise/damping.
    % Rotate phases so the maximum response floor becomes purely REAL.
    [~, max_floor] = max(abs(u_complex));
    phase_shift = angle(u_complex(max_floor));
    u_real = real(u_complex * exp(-1i * phase_shift));
    
    % Normalize mode shape relative to the top floor (Floor 6 = 1.0)
    phi_exp(:, m) = u_real / u_real(end);
end

% =========================================================================
% PLOT EXPERIMENTAL MODE SHAPES
% =========================================================================
floors = (1:num_channels)';
ground_floors = [0; floors]; % Include ground level (0) for plotting

figure('Name', 'FDD Experimental Mode Shapes', 'Color', 'w');

for m = 1:length(natural_freqs)
    subplot(1, length(natural_freqs), m);
    
    % Include ground (displacement = 0) at floor level 0
    mode_plot = [0; phi_exp(:, m)];
    
    plot(mode_plot, ground_floors, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    xline(0, 'k--', 'LineWidth', 0.8); % Centerline axis
    
    grid on;
    xlabel('Normalized Displacement');
    ylabel('Floor Level');
    title(sprintf('Mode %d: %.2f Hz', m, natural_freqs(m)));
    yticks(0:num_channels);
    ylim([0, num_channels + 0.5]);
end
%% FINDING STIFFNESS FROM THAT 
%  CONVERT NATURAL FREQUENCIES TO RAD/S
omega = 2 * pi * natural_freqs; % (rad/s)
num_modes_found = length(natural_freqs);

% 3. CONSTRUCT OVERDETERMINED SYSTEM A*k = b
A_total = [];
b_total = [];

for r = 1:num_modes_found
    phi = phi_exp(:, r);          % Mode shape vector (6x1)
    w2 = omega(r)^2;               % Squared natural frequency (rad^2/s^2)
    
    % Force vector for mode r
    b_r = w2 * M * phi;
    
    % Coefficient matrix for story stiffnesses k = [k1; k2; k3; k4; k5; k6]
    A_r = zeros(num_channels, num_channels);
    
    % Floor 1 equation
    A_r(1, 1) = phi(1);
    A_r(1, 2) = phi(1) - phi(2);
    
    % Intermediate Floors 2 to 5
    for i = 2:(num_channels - 1)
        A_r(i, i)   = phi(i) - phi(i-1);
        A_r(i, i+1) = phi(i) - phi(i+1);
    end
    
    % Top Floor 6
    A_r(num_channels, num_channels) = phi(num_channels) - phi(num_channels-1);
    
    % Stack matrices
    A_total = [A_total; A_r];
    b_total = [b_total; b_r];
end
options = optimoptions('lsqlin', 'Display', 'off');
k_identified = lsqlin(A_total, b_total, [], [], [], [], zeros(num_channels, 1), [], [], options);

% =========================================================================
%  DISPLAY RESULTS & RECONSTRUCT STIFFNESS MATRIX
% =========================================================================
fprintf('\n==========================================\n');
fprintf('    IDENTIFIED STORY STIFFNESSES (N/m)    \n');
fprintf('==========================================\n');
for i = 1:num_channels
    fprintf('  Story %d Stiffness (k_%d):  %10.2f N/m\n', i, i, k_identified(i));
end
fprintf('==========================================\n\n');

% Reconstruct the 6x6 Stiffness Matrix K
K_identified = zeros(num_channels);
K_identified(1, 1) = k_identified(1) + k_identified(2);
K_identified(1, 2) = -k_identified(2);

for i = 2:(num_channels - 1)
    K_identified(i, i-1) = -k_identified(i);
    K_identified(i, i)   = k_identified(i) + k_identified(i+1);
    K_identified(i, i+1) = -k_identified(i+1);
end
K_identified(num_channels, num_channels-1) = -k_identified(num_channels);
K_identified(num_channels, num_channels)   = k_identified(num_channels);
% =========================================================================
%% Professional Plotting of Acceleration Data
% =========================================================================
figure('Name', 'Raw Experimental Acceleration Time-History', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Create a Tiled Layout for clarity (6 subplots, one for each floor)
tiledlayout(2, 3, 'TileSpacing', 'compact');

% Color Palette for professional look
colors = lines(6); 

for i = 1:6
    nexttile;
    plot(raw_time, raw_accel(:, i), 'Color', colors(i,:), 'LineWidth', 1);
    
    % Professional Labeling
    ylabel(['Floor ', num2str(i), ' (m/s^2)'], 'FontSize', 10, 'FontWeight', 'bold');
    grid on;
    set(gca, 'Box', 'on', 'FontSize', 9);
    
    % Only show X-label on the bottom plot
    if i > 3
        xlabel('Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
    else
        set(gca, 'XTickLabel', []);
    end
end

% Add a main title
sgtitle('Acceleration Response of 6-Story Shear Building (Experimental)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% =========================================================================
% Professional Plotting of Calculated Displacement Data
% =========================================================================
figure('Name', 'Calculated Displacement Time-History', 'Color', 'w', 'Position', [150, 150, 800, 600]);

% Create a Tiled Layout for 6 floors
tiledlayout(2, 3, 'TileSpacing', 'compact');

% Use the same color palette for visual consistency across figures
colors = lines(6); 

for i = 1:6
    nexttile;
    % Plotting the displacement column for floor i
    plot(raw_time, displacement(:, i), 'Color', colors(i,:), 'LineWidth', 1);
    
    % Labeling with displacement units (meters)
    ylabel(['Floor ', num2str(i), ' (m)'], 'FontSize', 10, 'FontWeight', 'bold');
    grid on;
    set(gca, 'Box', 'on', 'FontSize', 9);
    
    % Apply X-label only to the bottom row of plots (tiles 4, 5, and 6)
    if i >= 4
        xlabel('Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
    else
        set(gca, 'XTickLabel', []);
    end
end

% Add a main title for the displacement response
sgtitle('Displacement Response of 6-Story Shear Building (Integrated)', ...
        'FontSize', 14, 'FontWeight', 'bold');


%% cleaning of the data 
data_detrend = detrend(exp_accel);

% 2. Setup Low-Pass Filter Parameters
dt = 0.00024;        % Your sampling interval
fc = 100;             % Cut-off frequency (50 Hz)
RC = 1 / (2 * pi * fc);
alpha = dt / (dt + RC);

% 3. Apply the Filter (Forward and Backward for Zero-Phase)
% We run it forward, then backward to ensure the peaks don't shift in time.
data_clean = zeros(size(data_detrend));

for i = 1:6
    % Forward Pass
    temp = zeros(size(data_detrend,1), 1);
    temp(1) = data_detrend(1, i);
    for n = 2:length(data_detrend)
        temp(n) = alpha * data_detrend(n, i) + (1 - alpha) * temp(n-1);
    end
    
    % Backward Pass (to eliminate time-lag/phase shift)
    data_clean(end, i) = temp(end);
    for n = length(data_detrend)-1:-1:1
        data_clean(n, i) = alpha * temp(n) + (1 - alpha) * data_clean(n+1, i);
    end
end

% 4. Verification Plot
figure('Name', 'Cleaned Signal Comparison', 'Color', 'w');
plot(exp_time, data_detrend(:,6), 'r', 'LineWidth', 0.5); hold on;
plot(exp_time, data_clean(:,6), 'k', 'LineWidth', 1.5);
title('Top Floor Acceleration: Raw vs Cleaned');
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
legend('Raw (Noisy)', 'Cleaned (Filtered)');
grid on;


%% ASSUMING I ONLY GET THE PAST AND PRESENT ACCELRATION THEN 
raw_accel_ = raw_accel;
n_steps = length(raw_time);
n_dof   = size(raw_accel_,2);

vel  = zeros(n_steps,n_dof);
disp = zeros(n_steps,n_dof);

vel(1,:)  = 0;
disp(1,:) = 0;
bias = 0; bias_vel=0;
alpha = 0.005; alpha_vel=0.0015;

for i=2:n_steps

    dt = raw_time(i)-raw_time(i-1);
    bias = (1-alpha)*bias + alpha*raw_accel_(i,:);

    raw_accel_(i,:) = raw_accel_(i,:) - bias;
    a0 = raw_accel_(i-1,:);
    a1 = raw_accel_(i,:);

    % Velocity
    vel(i,:) = vel(i-1,:) + 0.5*dt*(a0+a1);

    % bias_vel = (1-alpha_vel)*bias_vel + alpha_vel*vel(i,:);
    % 
    % vel(i,:) = vel(i,:) - bias_vel;
    % Displacement
    disp(i,:) = disp(i-1,:) ...
              + dt*vel(i-1,:) ...
              + 0.25*dt^2*(a0+a1);

end

figure(500);
plot(raw_time,displacement(:,6),raw_time,disp(:,6));
figure(501);
plot(raw_time,velocity(:,6),raw_time,vel(:,6));
% Trim the data to start at the moment of the snap-back

% displacement=disp;
% velocity = vel;
%% We look at Signal 6 (usually the top floor) to find the release
[~, idx] = max(abs(displacement(:, 6))); 
%idx=1;    % for considering the whole accn without any cut 
exp_time = raw_time(idx:end) - raw_time(idx);
exp_disp = displacement(idx:end, :);
exp_vel = velocity(idx:end, :);
exp_accel = raw_accel(idx:end, :);
figure(501);
plot(exp_time,exp_accel(:,6));
figure(502);
plot(exp_time,exp_disp(:,6));
figure(503);
plot(exp_time,exp_vel(:,6));
%% GOLAY FILTER 
% Coefficients for a 3rd order SG filter
sg_coeffs = [-0.0909; 0.0606; 0.1688; 0.2338; 0.2554; 0.2338; 0.1688; 0.0606; -0.0909];

% 2. Padding: This prevents the 'rounding' at the very start of the graph
% We pad with the first value so the filter sees a stable start
pad_n = 400000;
data_clean = zeros(size(data_detrend));

for i = 1:6
    % Add padding to the start and end of the signal
    signal_padded = [repmat(data_detrend(1,i), pad_n, 1); data_detrend(:,i); repmat(data_detrend(end,i), pad_n, 1)];
    
    % Apply filter using convolution
    filtered_padded = conv(signal_padded, sg_coeffs, 'same');
    
    % Strip the padding
    data_clean(:, i) = filtered_padded(pad_n+1 : end-pad_n);
end

% 3. Visualization
figure('Name', 'Cleaned Signal Comparison', 'Color', 'w');
plot(exp_time, data_detrend(:,6), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); hold on;
plot(exp_time, data_clean(:,6), 'k', 'LineWidth', 1.5);
title('Top Floor Acceleration: Savitzky-Golay Filtered');
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
legend('Raw (Noisy)', 'Cleaned (Preserved Peaks)');
grid on;

%% =========================================================================
% DAMPING CALCULATION (Logarithmic Decrement)
% =========================================================================
start_window=round(1.5/dt); end_window=round(5/dt);
data_clean=exp_disp(start_window:end_window,6);
% 1. Find the peaks in the top floor acceleration
min_height = 1.6e-5;       % Your baseline threshold

% Estimate the number of samples between two genuine consecutive crests 
% (e.g., if your wave period is roughly 50 samples, set this to ~35-40)
min_distance = 40;         

% Ensure the crest stands out from its local valleys to ignore shoulder ripples
min_prominence = 0.2e-5;   

% 3. Execute findpeaks with the optimized parameters
[peaks, locs] = findpeaks(data_clean, ...
    'MinPeakHeight', min_height, ...
    'MinPeakDistance', min_distance, ...
    'MinPeakProminence', min_prominence);

% 2. Select the first two consecutive peaks


% 3. Calculate Damping Ratio
delta = 0;
% Note: This assumes n=1 (consecutive peaks)
for i=2:size(peaks,1)
    x1 = peaks(i-1);
    x2 = peaks(i);
    delta = delta + log(x1 / x2);
   
end
 delta = delta/(size(peaks,1)-1); 
 zeta = delta / sqrt(4 * pi^2 + delta^2);
fprintf('Calculated Damping Ratio (zeta): %.4f\n', zeta);

% 4. Visualization of the Decay Envelope
figure('Color', 'w');
plot(exp_time(1:end_window-start_window+1), data_clean, 'k', 'LineWidth', 1); hold on;
plot(exp_time(locs), peaks, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
title(['Estimated Damping Ratio: \zeta = ', num2str(zeta, '%.4f')]);
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
grid on;



%% PROBLEM STATEMENT AND NUMERICAL MODEL
% True system parameters
m1 = 8.05; m2 = 8.06; m3 = 8.04; m4 = 7.95; m5 = 8.04; m6= 7.80;      % Masses (kg)
% As all of the columns have the same structure 
ko = 200000;
K_array_initial = [ 0.7706; 6.1796e+05 ; 5.1234e+05 ; 6.6990e+05;  6.1857e+05; 1.8292e+05 ];  % Initial spring stiffness (N/m)
alpha =2.1e-4; beta=8.6507e-04;       % Damping coefficients (Ns/m)
% for the we have to see rayleigh damping kibdof thing

% Mass matrix
M = diag([m1, m2, m3, m4, m5, m6]);

% Initial stiffness matrix (5x5 tridiagonal)
K_initial = build_stiffness_matrix(K_array_initial);
% Initial damping matrix (proportional damping, tridiagonal)
C_initial = alpha*M+beta*K_initial;

fprintf('Initial parameters:\n');
for l=1:6
    fprintf('initial stiffness k = %.3f \n', ...
        K_array_initial(l));
end
   fprintf('Rayleigh Damping Coefficient alpha = %.4f Ns/m, beta = %.4f Ns/m \n', ...
        alpha, beta);
% zonal stiffness (here its for spring )
Z = cell(6,1);

% Storey 1
Z{1} = zeros(6);
Z{1}(1,1) = 1;

% Storeys 2 to 6
for i = 2:6

    Zi = zeros(6);

    Zi(i-1,i-1) = 1;
    Zi(i-1,i)   = -1;
    Zi(i,i-1)   = -1;
    Zi(i,i)     = 1;

    Z{i} = Zi;

end
%% ---------------------------------------------------------
T = 5;              % Total simulation time (s)
N = round(T/dt);     % Number of time steps
time = linspace(0, T, N+1);

%% ---------------------------------------------------------
% Measurement Vectors ( size according to our iteration )***
% -------------------------------------------------------
accel_measured=exp_accel(1:20834,:)';
disp_measured=exp_disp(1:20834,:)';
vel_measured=exp_vel(1:20834,:)';
%% ---------------------------------------------------------
% DUAL ESTIMATION PARAMETERS 
% ----------------------------------------------------------
no_states = 12; % 6 disp, 6 vel
no_params = 8; % 6 damage , 2 damping coeff
no_obs_params = 6;  % 6 acceleraion measurements
no_obs_states = 2*no_obs_params; % disp and vel integrated from it 

% Storing History 
state_history = zeros(no_states,N+1);
param_history = zeros(no_params,N+1);

% Covarience matrix 
P_state  = 1.951605e-01*blkdiag(1,1,1,1,1,1,1,1,1,1,1,1); % 12x12 diag matrix 
P_param =  1.126978e-01*blkdiag(1,1,1,1,1,1,  1, 1);
%  q-5, R 1
% HAve to tune this part
Q_state = 1.000000e-18*eye(no_states);
Q_param =  4.233282e-14*eye(no_params);
R_state = 5.908816e-2*eye(no_obs_states);
R_param = 1.000000e+4*eye(no_obs_params);
r =5.228416e-01;

H_state = eye(no_states);  % remain same for entire iteraion
H_fd = zeros(no_params-2); % 6 fd and 2 damping coeff  
% Initialize the states and params
fd = zeros(6,1); taus = [alpha;beta]; x_param = [fd ;taus];
x_state = [0; 0; 0; 0.00001; 0.00008; 0.0001; 0; 0; 0; 0; 0; 0];
state_history(:,1)=x_state;

% Starting the DKF LOOP

for i=2:N+1
    % Building Kest for damage factor;
    for k = 1:6
        K_array(k) = K_array_initial(k)* (1 - fd(k));
    end
    K=build_stiffness_matrix(K_array);
    C=taus(1)*M + taus(2)*K;  % or current K???;
    
    %------- State Filter -----  

    A_cont = [zeros(6), eye(6); -M\K, -M\C];
    %B_cont = [zeros(6); inv(M)];
    A_disc = expm(A_cont * dt);
    %B_disc = (A_cont\(A_disc - eye(no_states))) * B_cont;

    x_state = A_disc*x_state;
    P_state = A_disc * P_state * A_disc'  + Q_state;
    
    Z_state = [disp_measured(:,i);vel_measured(:,i)];
    innovation_state = Z_state - H_state * x_state;
    S_state = H_state * P_state * H_state' + R_state;
    K_state = P_state * H_state' / S_state;
    x_state = x_state + K_state * innovation_state;
    P_state = (eye(no_states) - K_state * H_state) * P_state;  
    
    %------- Parameter Filter ---------------
    x_param=x_param + 1e-7*ones(8,1); %random walk update 
    P_param = P_param + Q_param;
    acc = accel_measured(:,i); %(x_state(7:12)-state_history(7:12,i-1))/dt; % or measured acc(:,i)
    G = M + dt*(1-r)*C + dt^2*(1-r)^2*K;
    % Observation matrix calculation for fds
    for k=1:6
        H_fd(:,k)= -K_array_initial(k).*Z{k}*(x_state(1:6)+dt*x_state(7:12)+dt^2*r*(1-r)*acc);
    end
    % Observation for alpha and beta
    H7 = M*(x_state(7:12)+dt*r*acc); 
    H8 = K*(x_state(7:12)+dt*r*acc);

    H_param = -G\[H_fd,H7,H8]; %6x8

    acc_ = -G\(K*(x_state(1:6)+dt*x_state(7:12)+dt^2*r*(1-r)*acc) + C*(x_state(7:12)+dt*r*acc));
    
    innovation_param= accel_measured(:,i) - acc_;
    S_param_ = H_param * P_param * H_param' + R_param;
    K_param_ = P_param * H_param' / S_param_;
    x_param = x_param + K_param_ * innovation_param;
    P_param = (eye(no_params) - K_param_ * H_param) * P_param; 

    fd=x_param(1:6); taus=x_param(7:8);
    % Constrain damage factors to [0, 0.8]
    fd = max(fd, 0);
    fd = min(fd, 0.99);
    % Constrain damage factors to [0, 0.2]
    taus(1) = max(taus(1), 0.00005);
    taus(1) = min(taus(1), 0.0005);
    taus(2) = max(taus(2), 0.000001);
    taus(2) = min(taus(2), 0.0001);

    % taus = max(taus, 0.001);
    % taus = min(taus, 0.2);
    x_param = [fd; taus];

    % Storing the variable 
    state_history(:,i)=x_state;
    param_history(:,i)=x_param;
end

%% plotting the damage factors 

figure('Name', 'Damage Factor Time-History', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Create a Tiled Layout for clarity (6 subplots, one for each floor)
tiledlayout(2, 3, 'TileSpacing', 'compact');

% Color Palette for professional look
colors = lines(6); 

for i = 1:6
    nexttile;
    plot(time, param_history(i,:), 'Color', colors(i,:), 'LineWidth', 1);
    
    % Professional Labeling
    ylabel(['Floor ', num2str(i), ' (unit)'], 'FontSize', 10, 'FontWeight', 'bold');
    grid on;
    set(gca, 'Box', 'on', 'FontSize', 9);
    
    % Only show X-label on the bottom plot
    if i > 3
        xlabel('Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
    else
        set(gca, 'XTickLabel', []);
    end
end

% Add a main title
sgtitle('Damage Response of 6-Story Shear Building (Experimental)', ...
        'FontSize', 14, 'FontWeight', 'bold');
%% pLOT THE DAMAGE FACTOR
figure(701);
plot(time,param_history(7,:),time,param_history(8,:));

figure(7002);
plot(time,state_history(6,:),exp_time,disp_measured(6,:));
%% ========================================================================
% HELPER FUNCTIONS
% ========================================================================
function K = build_stiffness_matrix(k)
    % Build stiffness matrix for 5-DOF system
    K = zeros(6,6);
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
    K(5,5) = k(5)+ k(6);
    K(5,6) = -k(6);
    K(6,5) = -k(6);
    K(6,6) = k(6);
end

%% uSING pso
% =========================================================================
% PSO TUNING FOR DKF HYPERPARAMETERS
% =========================================================================

% 1. Define bounds for the search space in log10 scale
% [log10(Q_state), log10(Q_param), log10(R_state), log10(R_param)]
% Example bounds: Q's between 1e-8 and 1e-1, R's between 1e-2 and 1e4
lb = [-4,   -4,   -8,  -8,  -8, -8,  0.2]; 
ub = [ 2,    2,    4,   4,    4,   4,    0.8];
 
% 2. Pack all required workspace variables into a structure 
% (Ensure these variables exist in your workspace before running)
sys_data.M = M;
sys_data.K_initial = K_initial;
sys_data.K_array_initial = K_array_initial;
sys_data.Z = Z;
sys_data.dt = dt;
sys_data.N = N;
sys_data.disp_measured = disp_measured;
sys_data.vel_measured = vel_measured;
sys_data.accel_measured = accel_measured;
sys_data.alpha = alpha;
sys_data.beta = beta;

% 3. Set up PSO options
options = optimoptions('particleswarm', ...
    'SwarmSize', 30, ...
    'MaxIterations', 20, ...
    'Display', 'iter', ...
    'UseParallel', false); % Set to true if you have Parallel Computing Toolbox

% 4. Run the PSO Optimization
objective_func = @(tune_vars) dkf_frequency_objective(tune_vars, sys_data);
optimal_logs = particleswarm(objective_func, 7, lb, ub, options);

% 5. Extract and display optimal linear values
opt_P_state = 10^optimal_logs(1);
opt_P_param = 10^optimal_logs(2);
opt_Q_state = 10^optimal_logs(3);
opt_Q_param = 10^optimal_logs(4);
opt_R_state = 10^optimal_logs(5);
opt_R_param = 10^optimal_logs(6);
opt_r=optimal_logs(7);
fprintf('\n--- OPTIMIZATION COMPLETE ---\n');
fprintf('Optimal P_state scalar: %e\n', opt_P_state);
fprintf('Optimal P_param scalar: %e\n', opt_P_param);
fprintf('Optimal Q_state scalar: %e\n', opt_Q_state);
fprintf('Optimal Q_param scalar: %e\n', opt_Q_param);
fprintf('Optimal R_state scalar: %e\n', opt_R_state);
fprintf('Optimal R_param scalar: %e\n', opt_R_param);
fprintf('Optimal r scalar: %e\n',opt_r);

function cost = dkf_frequency_objective(tune_vars, sys)
    % ---------------------------------------------------------
    % UNPACK SYSTEM VARIABLES
    % ---------------------------------------------------------
    M = sys.M; K_array_initial = sys.K_array_initial;Z = sys.Z; %ko = sys.ko;
    %  
    dt = sys.dt; N = sys.N; alpha = sys.alpha; beta = sys.beta;
    disp_measured = sys.disp_measured; vel_measured = sys.vel_measured;
    accel_measured = sys.accel_measured;
    
    % ---------------------------------------------------------
    % INITIALIZE DKF PARAMETERS
    % ---------------------------------------------------------
    no_states = 12; 
    no_params = 8; 
    no_obs_params = 6;  
    no_obs_states = 2*no_obs_params; 
    
    state_history = zeros(no_states, N+1);
    param_history = zeros(no_params, N+1);
    
    % P_state  = 1e-3 * blkdiag(1,1,1,1,1,1,1,1,1,1,1,1); 
    % P_param = 1e-2 * blkdiag(1,1,1,1,1,1,1,1);
    
    % Apply PSO trial parameters (converting from log10 back to linear)
    P_state = (10^tune_vars(1)) * eye(no_states);
    P_param = (10^tune_vars(2)) * eye(no_params);
    Q_state = (10^tune_vars(3)) * eye(no_states);
    Q_param = (10^tune_vars(4)) * eye(no_params);
    R_state = (10^tune_vars(5)) * eye(no_obs_states);
    R_param = (10^tune_vars(6)) * eye(no_obs_params);
    r = tune_vars(7);
   
    H_state = eye(no_states);  
    H_fd = zeros(no_params-2); % 6x6 matrix 
    
    fd = zeros(6,1); 
    taus = [alpha; beta]; 
    x_param = [fd; taus];
    x_state = [0; 0; 0; 0.00001; 0.00008; 0.0001; 0; 0; 0; 0; 0; 0];
    state_history(:,1) = x_state;
    K_array = zeros(1, 6);
    
    % ---------------------------------------------------------
    % START DKF LOOP
    % ---------------------------------------------------------
    for i = 2:N+1
        % Build K for damage factor
        for k = 1:6
            K_array(k) = K_array_initial(k) * (1 - fd(k));
        end
        % Assuming build_stiffness_matrix is on your path
        K = build_stiffness_matrix(K_array);
        C = taus(1)*M + taus(2)*K;  % Corrected to use current K
        
        % ------- State Filter -----  
        A_cont = [zeros(6), eye(6); -M\K, -M\C];
        A_disc = expm(A_cont * dt);
        
        x_state = A_disc * x_state;
        P_state = A_disc * P_state * A_disc' + Q_state;
        
        Z_state = [disp_measured(:,i); vel_measured(:,i)];
        innovation_state = Z_state - H_state * x_state;
        S_state = H_state * P_state * H_state' + R_state;
        K_state = P_state * H_state' / S_state;
        x_state = x_state + K_state * innovation_state;
        P_state = (eye(no_states) - K_state * H_state) * P_state;  
        
        % ------- Parameter Filter ---------------
        x_param = x_param + 1e-7*ones(8,1); 
        P_param = P_param + Q_param;
        acc = accel_measured(:,i); 
        G = M + dt*(1-r)*C + dt^2*(1-r)^2*K;
        
        for k = 1:6
            H_fd(:,k) = -K_array_initial(k) .* Z{k} * (x_state(1:6) + dt*x_state(7:12) + dt^2*r*(1-r)*acc);
        end
        
        H7 = M * (x_state(7:12) + dt*r*acc); 
        H8 = K * (x_state(7:12) + dt*r*acc);
        H_param = -G \ [H_fd, H7, H8]; 
        
        acc_ = -G \ (K * (x_state(1:6) + dt*x_state(7:12) + dt^2*r*(1-r)*acc) + C * (x_state(7:12) + dt*r*acc));
        
        innovation_param = accel_measured(:,i) - acc_;
        S_param_ = H_param * P_param * H_param' + R_param;
        K_param_ = P_param * H_param' / S_param_;
        x_param = x_param + K_param_ * innovation_param;
        P_param = (eye(no_params) - K_param_ * H_param) * P_param; 
        
        fd = x_param(1:6); 
        taus = x_param(7:8);
        
        % Constraints
        fd = max(fd, -4);
        fd = min(fd, 0.99);
        taus(1) = max(taus(1), 0.00005);
        taus(1) = min(taus(1), 0.0001);
        taus(2) = max(taus(2), 0.000001);
        taus(2) = min(taus(2), 0.00001);
        x_param = [fd; taus];
        
        % Store
        state_history(:,i) = x_state;
        param_history(:,i) = x_param;
    end
    
    % ---------------------------------------------------------
    % EVALUATE FINAL FREQUENCIES (COST FUNCTION)
    % ---------------------------------------------------------
    % Reconstruct the final stiffness matrix
    fd_final = param_history(1:6, end);
    K_array_final = zeros(1, 6);
    for k = 1:6
        K_array_final(k) = K_array_initial(k) * (1 - fd_final(k));
    end
    K_final = build_stiffness_matrix(K_array_final);
    
    % Solve generalized eigenvalue problem: K*v = lambda*M*v
    [~, D] = eig(K_final, M);
    
    % Extract eigenvalues (lambda = omega^2)
    eigenvalues = diag(D);
    
    % Ensure no negative/complex eigenvalues due to severe constraint violations
    if any(eigenvalues <= 0) || any(~isreal(eigenvalues))
        cost = 1e6; % Heavy penalty if model becomes non-physical
        return;
    end
    
    % Convert to Hertz: f = sqrt(lambda) / (2*pi)
    frequencies_Hz = sqrt(eigenvalues) / (2*pi);
    frequencies_Hz = sort(frequencies_Hz); % Ascending order
    
    % Calculate RMSE against target frequencies
    target_freqs = [6.1; 20.31; 34.55; 55.0; 71.18; 83.4];
    calculated_freqs = frequencies_Hz(1:6);
    
    % Cost is the norm (or RMSE) of the difference
    cost = norm(calculated_freqs - target_freqs); 
end


%% pso FOR EKF 
% 1. Define bounds for the search space in log10 scale
% Optimization variables: [log10(P_aug), log10(Q_aug), log10(R_meas)]
% We optimize scalars multiplying identity matrices for the augmented system
lb = [-8,   -8,   -8.5]; 
ub = [ 2,    2,    2.0];
 
% 2. Pack all required workspace variables into a structure 
sys_data.M = M;
sys_data.K_array_initial = K_array_initial;
sys_data.ko = ko;
sys_data.Z = Z;
sys_data.dt = dt;
sys_data.N = N;
sys_data.disp_measured = disp_measured;
sys_data.vel_measured = vel_measured;
sys_data.accel_measured = accel_measured;
sys_data.alpha = alpha;
sys_data.beta = beta;

% 3. Set up PSO options
options = optimoptions('particleswarm', ...
    'SwarmSize', 30, ...
    'MaxIterations', 20, ...
    'Display', 'iter', ...
    'UseParallel', false); 

% 4. Run the PSO Optimization (Optimizing 3 tuning variables)
objective_func = @(tune_vars) ekf_frequency_objective(tune_vars, sys_data);
optimal_logs = particleswarm(objective_func, 3, lb, ub, options);

% 5. Extract and display optimal linear values
opt_P_aug = 10^optimal_logs(1);
opt_Q_aug = 10^optimal_logs(2);
opt_R_meas = 10^optimal_logs(3);

fprintf('\n--- EKF OPTIMIZATION COMPLETE ---\n');
fprintf('Optimal Augmented P scalar: %e\n', opt_P_aug);
fprintf('Optimal Augmented Q scalar: %e\n', opt_Q_aug);
fprintf('Optimal Measurement R scalar: %e\n', opt_R_meas);

function cost = ekf_frequency_objective(tune_vars, sys)
    % ---------------------------------------------------------
    % UNPACK SYSTEM VARIABLES
    % ---------------------------------------------------------
    M = sys.M; K_array_initial = sys.K_array_initial; ko = sys.ko; Z = sys.Z;
    dt = sys.dt; N = sys.N; alpha = sys.alpha; beta = sys.beta;
    disp_measured = sys.disp_measured; vel_measured = sys.vel_measured;
    
    % ---------------------------------------------------------
    % INITIALIZE AUGMENTED EKF PARAMETERS
    % ---------------------------------------------------------
    no_states = 12;  % Displacement (6) + Velocity (6)
    no_params = 6;   % Structural damage parameters (fd)
    no_aug = no_states + no_params; % Total augmented state dimension (18)
    
    no_obs_states = 12; % Measured displacements (6) + velocities (6)
    
    aug_history = zeros(no_aug, N+1);
    
    % Apply PSO trial parameters (converting from log10 back to linear)
    P_aug = (10^tune_vars(1)) * eye(no_aug);
    Q_aug = (10^tune_vars(2)) * eye(no_aug);
    R_meas = (10^tune_vars(3)) * eye(no_obs_states);
    
    % Initial states initialization
    x_state = [0; 0; 0; 0.00001; 0.00008; 0.0001; 0; 0; 0; 0; 0; 0];
    fd = zeros(6, 1); 
    
    % Formulate initial augmented state vector: [disp; vel; fd]
    x_aug = [x_state; fd];
    aug_history(:, 1) = x_aug;
    
    K_array = zeros(1, 6);
    invM = inv(M);
    
    % ---------------------------------------------------------
    % START EKF STATE-ESTIMATION LOOP
    % ---------------------------------------------------------
    for i = 2:N+1
        % Extract current estimates from augmented vector
        u    = x_aug(1:6);
        v    = x_aug(7:12);
        fd_c = x_aug(13:18);
        
        % 1. Reconstruct current parameter-dependent matrices
        for k = 1:6
            K_array(k) = K_array_initial(k) * (1 - fd_c(k));
        end
        K = build_stiffness_matrix(K_array);
        C = alpha*M + beta*K;
        
        % 2. State Transition Jacobian (A_matrix) calculation via continuous-to-discrete
        % Derivative of state dynamics with respect to fd parameters
        dF_dfd = zeros(6, 6);
        for k = 1:6
            % Z{k} maps structural layout to stiffness components
            dF_dfd = dF_dfd - K_array_initial(k) * Z{k} * (u + beta*v); 
        end
        
        % Continuous Jacobian evaluation
        Ac = [  zeros(6),             eye(6),         zeros(6,6);
              -invM * K,           -invM * C,       invM * dF_dfd;
                zeros(6,6),           zeros(6,6),     zeros(6,6)  ];
            
        A_disc = expm(Ac * dt);
        
        % 3. EKF Time Update (Prediction Step)
        % State equations for random-walk parameters imply parameter drift remains zero
        u_pred = u + dt*v;
        v_pred = v + dt * (invM * (-K*u - C*v));
        fd_pred = fd_c; 
        
        x_aug = [u_pred; v_pred; fd_pred];
        P_aug = A_disc * P_aug * A_disc' + Q_aug;
        
        % 4. EKF Measurement Update (Correction Step)
        H_meas = [eye(no_obs_states), zeros(no_obs_states, no_params)]; 
        
        Z_actual = [disp_measured(:,i); vel_measured(:,i)];
        innovation = Z_actual - H_meas * x_aug;
        
        S_matrix = H_meas * P_aug * H_meas' + R_meas;
        K_gain = P_aug * H_meas' / S_matrix;
        
        x_aug = x_aug + K_gain * innovation;
        P_aug = (eye(no_aug) - K_gain * H_meas) * P_aug;
        
        % 5. Enforce Physical Constraints Post-Correction
        x_aug(13:18) = max(x_aug(13:18), -1);
        x_aug(13:18) = min(x_aug(13:18), 0.99);
        
        % Store trace
        aug_history(:, i) = x_aug;
    end
    
    % ---------------------------------------------------------
    % EVALUATE FINAL FREQUENCIES (COST FUNCTION)
    % ---------------------------------------------------------
    fd_final = aug_history(13:18, end);
    K_array_final = zeros(1, 6);
    for k = 1:6
        K_array_final(k) = K_array_initial(k) * (1 - fd_final(k));
    end
    K_final = build_stiffness_matrix(K_array_final);
    
    % Solve generalized eigenvalue problem: K*v = lambda*M*v
    [~, D] = eig(K_final, M);
    eigenvalues = diag(D);
    
    % Heavy penalty if model equations become non-physical
    if any(eigenvalues <= 0) || any(~isreal(eigenvalues))
        cost = 1e6; 
        return;
    end
    
    % Convert to Hertz and calculate RMSE error metrics
    frequencies_Hz = sort(sqrt(eigenvalues) / (2*pi));
    target_freqs = [6.1; 20.3; 34.6; 55; 71.2; 83.4];
    calculated_freqs = frequencies_Hz(1:6);
    
    cost = norm(calculated_freqs - target_freqs); 
end

%% =========================================================================
% =========================================================================
% PSO TUNING FOR DKF HYPERPARAMETERS (DIRECT STIFFNESS TRACKING)
% =========================================================================

% 1. Define bounds for search space in log10 scale
lb = [-2,   9,   -5,  -23, 11,  4,  0.7]; 
ub = [  -1,  12,  -4, -21,  13, 6,  0.8];
 
% 2. Pack required workspace variables into structure
sys_data.M = M;
sys_data.K_initial = K_initial;
sys_data.K_array_initial = K_array_initial;
sys_data.Z = Z;
sys_data.dt = dt;
sys_data.N = N;
sys_data.disp_measured = disp_measured;
sys_data.vel_measured = vel_measured;
sys_data.accel_measured = accel_measured;
sys_data.alpha = 1.336;
sys_data.beta = 6.65e-5;

% 3. Set up PSO options
options = optimoptions('particleswarm', ...
    'SwarmSize', 4, ...
    'MaxIterations', 3, ...
    'Display', 'iter', ...
    'UseParallel', false); 

% 4. Run PSO Optimization
objective_func = @(tune_vars) dkf_frequency_objective_direct_stiff(tune_vars, sys_data);
optimal_logs = particleswarm(objective_func, 7, lb, ub, options);

% 5. Extract optimal values
opt_P_state = 10^optimal_logs(1);
opt_P_param = 10^optimal_logs(2);
opt_Q_state = 10^optimal_logs(3);
opt_Q_param = 10^optimal_logs(4);
opt_R_state = 10^optimal_logs(5);
opt_R_param = 10^optimal_logs(6);
opt_r       = optimal_logs(7);

fprintf('\n--- OPTIMIZATION COMPLETE FOR DIRECT STIFFNESS IN DKF ---\n');
fprintf('Optimal P_state scalar: %e\n', opt_P_state);
fprintf('Optimal P_param scalar: %e\n', opt_P_param);
fprintf('Optimal Q_state scalar: %e\n', opt_Q_state);
fprintf('Optimal Q_param scalar: %e\n', opt_Q_param);
fprintf('Optimal R_state scalar: %e\n', opt_R_state);
fprintf('Optimal R_param scalar: %e\n', opt_R_param);
fprintf('Optimal r scalar:       %e\n', opt_r);


% OBJECTIVE FUNCTION
function cost = dkf_frequency_objective_direct_stiff(tune_vars, sys)
    % Unpack System Data
    M = sys.M; K_array_initial = sys.K_array_initial; Z = sys.Z;
    dt = sys.dt; N = sys.N; alpha = 1.336; beta = 6.665e-5;
    disp_measured = sys.disp_measured; 
    vel_measured  = sys.vel_measured;
    accel_measured = sys.accel_measured;
    
    no_states = 12; 
    no_params = 8; 
    no_obs_params = 6;  
    no_obs_states = 2 * no_obs_params; 
    
    param_history = zeros(no_params, N+1);
    
    % Initialize Covariance Matrices Once
    P_state = (10^tune_vars(1)) * eye(no_states);
    P_param = (10^tune_vars(2)) * eye(no_params);
    Q_state = (10^tune_vars(3)) * eye(no_states);
    Q_param = (10^tune_vars(4)) * eye(no_params);
    R_state = (10^tune_vars(5)) * eye(no_obs_states);
    R_param = (10^tune_vars(6)) * eye(no_obs_params);
    r       = tune_vars(7);
   
    H_state = eye(no_states);  
    H_K = zeros(no_obs_params, 6);
    
    K_vals = K_array_initial; 
    taus = [alpha; beta]; 
    x_param = [K_vals; taus];
    x_state = zeros(no_states, 1);
    x_state(4:6) = [1e-5; 8e-5; 1e-4]; % Initial non-zero states
    
    I_state = eye(no_states);
    I_param = eye(no_params);
    
    K_array = zeros(1, 6);
    
    % ---------------------------------------------------------
    % DUAL KALMAN FILTER LOOP
    % ---------------------------------------------------------
    for i = 2:N+1
        K_array = x_param(1:6)';
        K = build_stiffness_matrix(K_array);
        C = x_param(7)*M + x_param(8)*K;  
        
        % --- State Filter ---
        A_cont = [zeros(6), eye(6); -M\K, -M\C];
        A_disc = expm(A_cont * dt);
        
        x_state = A_disc * x_state;
        P_state = A_disc * P_state * A_disc' + Q_state; % P_state updates recursively
        
        Z_state = [disp_measured(:,i-1); vel_measured(:,i-1)];
        innovation_state = Z_state - H_state * x_state;
        S_state = H_state * P_state * H_state' + R_state;
        K_state = P_state * H_state' / S_state;
        x_state = x_state + K_state * innovation_state;
        P_state = (I_state - K_state * H_state) * P_state;  
        
        % --- Parameter Filter ---
        % (Removed artificial 1e-7 drift addition)
        P_param = P_param + Q_param;                   % P_param updates recursively
        acc = accel_measured(:,i-1); 
        G = M + dt*(1-r)*C + dt^2*(1-r)^2*K;
        
        for k = 1:6
            H_K(:,k) = Z{k} * (x_state(1:6) + dt*x_state(7:12) + dt^2*r*(1-r)*acc);
        end
        
        H7 = M * (x_state(7:12) + dt*r*acc); 
        H8 = K * (x_state(7:12) + dt*r*acc);
        H_param = -G \ [H_K, H7, H8]; 
        
        acc_ = -G \ (K * (x_state(1:6) + dt*x_state(7:12) + dt^2*r*(1-r)*acc) + C * (x_state(7:12) + dt*r*acc));
        
        innovation_param = acc - acc_;
        S_param_ = H_param * P_param * H_param' + R_param;
        K_param_ = P_param * H_param' / S_param_;
        x_param = x_param + K_param_ * innovation_param;
        P_param = (I_param - K_param_ * H_param) * P_param; 
        
        % Physical Bounding
        x_param(1:6) = max(min(x_param(1:6), 5.0 * 200000), 0.01 * 20000);
        x_param(7:8) = max(min(x_param(7:8), 0.5), 0.0002);
        
        param_history(:,i) = x_param;
    end
    
    % ---------------------------------------------------------
    % COST EVALUATION (Averaged over last 20% of time history)
    % ---------------------------------------------------------
    stable_indices = round(0.8 * N):N+1;
    K_vals_mean = mean(param_history(1:6, stable_indices), 2);
    
    for l=1:6
        fprintf('K value: %e\n', K_vals_mean(l));
    end
    fprintf('alpha value: %e and beta value: %e  \n', param_history(7,end),param_history(8,end));

    K_final = build_stiffness_matrix(K_vals_mean');
    
    [~, D] = eig(K_final, M);
    eigenvalues = diag(D);
    
    if any(eigenvalues <= 0) || any(~isreal(eigenvalues))
        cost = 1e6; 
        return;
    end
    
    frequencies_Hz = sort(sqrt(eigenvalues) / (2*pi));
    
    % Target Frequencies (Matching 6-DOF dimension)
    target_freqs = [6.1; 20.3; 34.6; 55; 71.2; 83.4]; 
    
    % Objective: Frequency error + parameter variance penalty
    freq_error = norm(frequencies_Hz - target_freqs);
    param_var_penalty = sum(std(param_history(1:6, stable_indices), 0, 2) / 200000);
    
    cost = freq_error + 0.1 * param_var_penalty;
end
%% uSING MINCON FUNCTION OPTIMIZATIONN 
m = [8.05; 8.06; 8.04; 7.95; 8.04; 7.80]; % Storey masses (tons)
M = diag(m);

% Measured Targets
f_target = [6.1; 20.31; 34.55; 55.0; 71.18; 83.4]; % Target 5 frequencies (Hz)

% Load/Define Experimental Acceleration Data
% Assume: t_meas = [N x 1] time vector, acc_meas = [N x 6] matrix for 6 floors
% ground_acc = [N x 1] ground motion input (e.g., base excitation / impulse)
% (For demonstration, ensure these variables are loaded in your workspace)
% load('measured_vibration_data.mat'); 

% Damping matrix assumption (e.g., Rayleigh damping or 2% modal damping)
zeta = 0.02; % 2% damping ratio

% 2. Weighting Factors
w_freq = 1.0; 
w_acc  = 0.5;
acc_meas = accel_measured';
window_size = 100;
% 3. Hybrid Objective Function
objective = @(k) hybrid_objective(k, M, f_target, time, acc_meas, w_freq, w_acc,window_size);

% 4. Optimization Setup
k0 = 200000 * ones(6, 1);
lb = 1e3 * ones(6, 1);
ub = 9e8 * ones(6, 1);

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'FiniteDifferenceType', 'central', ...
    'FiniteDifferenceStepSize', 1e-3, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12, ...
    'MaxFunctionEvaluations', 2000000);

% 5. Run Optimization
tic;
[k_opt, fval, exitflag, output] = fmincon(objective, k0, [], [], [], [], lb, ub, [], options);
execution_time = toc;

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function J = hybrid_objective(k, M, f_target, t, acc_meas, w_f, w_a, window_size)
    % --- Part A: Frequency Error ---
    f_sim = get_freqs(k, M, length(f_target));
    err_freq = norm(f_sim - f_target);
    
    % --- Part B: Window-Averaged Acceleration Error ---
    acc_sim = simulate_accelerations(k, M, t); % Size: [N_steps x N_floors]
    
    num_steps = size(acc_sim, 1);
    num_floors = size(acc_sim, 2);
    
    % Ensure the time steps can be evenly divided into windows
    num_windows = floor(num_steps / window_size);
    valid_steps = num_windows * window_size;
    
    % Trim arrays to fit whole windows if total steps isn't a multiple of window_size
    acc_sim_trimmed = acc_sim(1:valid_steps, :);
    acc_meas_trimmed = acc_meas(1:valid_steps, :);
    
    err_acc = 0;
    
    for j = 1:num_floors
        % Reshape to [window_size x num_windows] and calculate column-wise means
        % Resulting vectors have length = num_windows (e.g., 20,000 / 100 = 200 points)
        sim_win_avg  = mean(reshape(acc_sim_trimmed(:, j), window_size, num_windows), 1)';
        meas_win_avg = mean(reshape(acc_meas_trimmed(:, j), window_size, num_windows), 1)';
        
        % Calculate L2-norm on the windowed averages
        err_acc = err_acc + norm(sim_win_avg - meas_win_avg) / norm(meas_win_avg);
    end
    
    err_acc = err_acc / num_floors;
    
    % Combined Objective
    J = w_f * err_freq + w_a * err_acc;
end

function f_firstFEQ = get_freqs(k, M, num_freqs)
    K = assemble_K(k);
    [~, D] = eig(K, M);
    eigenvalues = real(diag(D));
    eigenvalues(eigenvalues < 0) = 0;
    freqs = sort(sqrt(eigenvalues) / (2 * pi));
    f_firstFEQ = freqs(1:num_freqs);
end

function K = assemble_K(k)
    K = zeros(6, 6);
    K(1,1) = k(1) + k(2);
    K(1,2) = -k(2);
    for i = 2:5
        K(i, i-1) = -k(i);
        K(i, i)   = k(i) + k(i+1);
        K(i, i+1) = -k(i+1);
    end
    K(6,5) = -k(6);
    K(6,6) = k(6);
end

function acc_sim = simulate_accelerations(k, M, t, zeta)
% NEWMARK-BETA FREE VIBRATION ACCELERATION SIMULATION
% Solves M*x'' + C*x' + K*x = 0 using Average Acceleration Method (gamma=0.5, beta=0.25)

    N = length(k);            % 6 DOFs
    num_steps = length(t);
    dt = 0.00024;        % Sampling interval
    
    % Assemble Stiffness Matrix K
    K = assemble_K(k);
    
    % Fixed Rayleigh Damping Coefficients
    alpha = 2.154170e-04;
    beta_rayleigh = 9.154170e-04;
    C = alpha * M + beta_rayleigh * K;

    % Newmark Integration Parameters (Average Acceleration)
    gamma = 0.5;
    beta  = 0.25;

    % Pre-allocate State Vectors
    x     = zeros(N, num_steps); % Displacement
    x_dot = zeros(N, num_steps); % Velocity
    x_ddot = zeros(N, num_steps); % Acceleration
    
    % Initial Conditions at t = 0
    x(:, 1)     = 0.0000001; %u0;            % Initial floor displacements (m)
    x_dot(:, 1) = 0.0000001; %v0;            % Initial floor velocities (m/s)

    % Initial Acceleration at t = 0 from Equilibrium: M*x_ddot0 + C*x_dot0 + K*x0 = 0
    x_ddot(:, 1) = M \ (-C * x_dot(:, 1) - K * x(:, 1));

    % Integration Integration Constants
    a0 = 1 / (beta * dt^2);
    a1 = gamma / (beta * dt);
    a2 = 1 / (beta * dt);
    a3 = 1 / (2 * beta) - 1;
    a4 = (gamma / beta) - 1;
    a5 = (dt / 2) * ((gamma / beta) - 2);

    % Effective Stiffness Matrix (K_hat = K + a0*M + a1*C)
    K_hat = K + a0 * M + a1 * C;

    % Time-Stepping Loop (Free Vibration)
    for i = 1:(num_steps - 1)
        % Effective Load Vector at step i+1 (Dependence on past states only)
        P_eff = M * (a0 * x(:, i) + a2 * x_dot(:, i) + a3 * x_ddot(:, i)) + ...
                C * (a1 * x(:, i) + a4 * x_dot(:, i) + a5 * x_ddot(:, i));

        % Solve for Next Step Displacement
        x(:, i+1) = K_hat \ P_eff;

        % Update Accelerations and Velocities
        x_ddot(:, i+1) = a0 * (x(:, i+1) - x(:, i)) - a2 * x_dot(:, i) - a3 * x_ddot(:, i);
        x_dot(:, i+1)  = x_dot(:, i) + dt * ((1 - gamma) * x_ddot(:, i) + gamma * x_ddot(:, i+1));
    end

    % Return acceleration matrix (size: length(t) x N)
    acc_sim = x_ddot';
    %size(acc_sim)
end

%% Extract and Sort Mode Shapes and Frequencies
k_opt_matrix = build_stiffness_matrix(k_opt);
[PHIS, EIGS] = eig(k_opt_matrix, M);

% 1. Extract natural frequencies (rad/s and Hz)
omegas = sqrt(diag(EIGS));
[omegas_sorted, sort_idx] = sort(omegas); % Sort in ascending order
freqs_Hz = omegas_sorted / (2 * pi);

% 2. Rearrange mode shape matrix corresponding to sorted frequencies
phi_sorted = PHIS(:, sort_idx);

% 3. Include Ground Level (Floor 0 = 0 displacement) for plotting
num_floors = 6;
floors = 0:num_floors; % Floors 0 to 6

% Pre-allocate mode matrix including base (7 x 6)
phi_plot = zeros(num_floors + 1, num_floors);

for mode = 1:num_floors
    % Extract mode shape vector for current mode
    v = phi_sorted(:, mode);
    
    % Normalize by top floor displacement (so top floor is always +1)
    v_norm = v / v(end);
    
    % Store with zero displacement at floor 0
    phi_plot(:, mode) = [0; v_norm];
end

% Plot First 5 Mode Shapes
figure('Name', 'Mode Shapes of 6-Story Building', 'Color', 'w');

colors = ['r', 'b', 'g','m','k'];
markers = ['o', 's', '^','*','.'];

for mode = 1:5
    subplot(1, 5, mode);
    
    % Draw vertical reference line (zero displacement)
    plot([0 0], [0 num_floors], 'k--', 'LineWidth', 1); hold on;
    
    % Plot mode shape profile
    plot(phi_plot(:, mode), floors, ...
        'Color', colors(mode), ...
        'LineWidth', 2, ...
        'Marker', markers(mode), ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', colors(mode));
    
    % Formatting
    grid on;
    title(sprintf('Mode %d (%.2f Hz)', mode, freqs_Hz(mode)), 'FontSize', 11);
    xlabel('Normalized Displacement');
    ylabel('Floor Level');
    ylim([0 num_floors]);
    yticks(0:num_floors);
    yticklabels({'Base (0)', 'Floor 1', 'Floor 2', 'Floor 3', 'Floor 4', 'Floor 5', 'Floor 6'});
    xlim([-1.5, 1.5]);
end

sgtitle('First 6 Vibration Mode Shapes (fmincon Optimized k)', 'FontSize', 14, 'FontWeight', 'bold');
%% =========================================================================
% NEWMARK-BETA FREE VIBRATION ACCELERATION HISTORY
% =========================================================================

% 1. System Parameters
                % Time step (seconds)
t_final = 5.0;               % Simulation duration (seconds)
time = 0:dt:t_final;
N_pts = length(time);

% Rayleigh Damping Coefficients (e.g., 2% damping ratio)
alpha_rayleigh = 2.000000e-04; 
beta_rayleigh  = 9.154170e-04 ;

% 2. Build System Matrices using optimized stiffness (k_opt)
M_mat = diag(m);
K_mat = build_stiffness_matrix([8.063967e+04; 6.026747e+05 ; 5.531480e+05 ; 6.425025e+05;6.258761e+05;  1.770501e+05]);
C_mat = alpha_rayleigh * M_mat + beta_rayleigh * K_mat;

% 3. Initial Conditions (10 cm initial displacement at top floor)
u0 = zeros(6, 1);
u0(6) = 0.00652;                % 0.1 m displacement on Floor 6
v0 = zeros(6, 1);            % Initial velocity = 0

% Compute Initial Acceleration a0 from M*a0 + C*v0 + K*u0 = 0
a0 = M_mat \ (-C_mat * v0 - K_mat * u0);

% 4. Newmark-Beta Integration Setup (\gamma = 0.5, \beta = 0.25)
gamma = 0.5; 
beta_nb = 0.25;

a0_c = 1 / (beta_nb * dt^2);
a1_c = gamma / (beta_nb * dt);
a2_c = 1 / (beta_nb * dt);
a3_c = 1 / (2 * beta_nb) - 1;
a4_c = gamma / beta_nb - 1;
a5_c = (dt / 2) * (gamma / beta_nb - 2);

% Effective Stiffness Matrix (constant across time steps)
K_eff = K_mat + a0_c * M_mat + a1_c * C_mat;

% Allocate Time History Matrices (6 x N_pts)
u_hist = zeros(6, N_pts);
v_hist = zeros(6, N_pts);
a_hist = zeros(6, N_pts);

% Set Step 1
u_hist(:, 1) = u0;
v_hist(:, 1) = v0;
a_hist(:, 1) = a0;

% 5. Time Integration Loop
for k = 1:N_pts-1
    u_k = u_hist(:, k);
    v_k = v_hist(:, k);
    a_k = a_hist(:, k);
    
    % Effective Force
    F_eff = M_mat * (a0_c * u_k + a2_c * v_k + a3_c * a_k) + ...
            C_mat * (a1_c * u_k + a4_c * v_k + a5_c * a_k);
        
    % Solve for u_{k+1}
    u_next = K_eff \ F_eff;
    
    % Update acceleration and velocity
    a_next = a0_c * (u_next - u_k) - a2_c * v_k - a3_c * a_k;
    v_next = v_k + a5_c * a_k + (gamma * dt) * a_next;
    
    % Store
    u_hist(:, k+1) = u_next;
    v_hist(:, k+1) = v_next;
    a_hist(:, k+1) = a_next;
end

% 6. Plot Free Vibration Acceleration Histories
figure('Name', 'Free Vibration Acceleration History', 'Color', 'w');
colors = lines(6);

for i = 1:6
    subplot(2, 3, i);
    plot(time, a_hist(i, :), 'Color', colors(i, :), 'LineWidth', 1.2);
    grid on;
    title(sprintf('Floor %d Acceleration', i));
    xlabel('Time (s)');
    ylabel('Accel (m/s^2)');
    xlim([0, t_final]);
end

% %% Helper Function to Construct Tri-Diagonal Stiffness Matrix
% function K = build_stiffness_matrix(k)
%     K = zeros(6,6);
%     K(1,1) = k(1) + k(2); 
%     K(1,2) = -k(2);
%     for j = 2:5
%         K(j, j-1) = -k(j);
%         K(j, j)   = k(j) + k(j+1);
%         K(j, j+1) = -k(j+1);
%     end
%     K(6,5) = -k(6); 
%     K(6,6) = k(6);
% end
%% ---------------------------------------------------------
% EIGENVALUE SOLUTION
% ----------------------------------------------------------
% for k = 1:6
%         K_array(k) = K_array_initial(k) * (1 - fd(k));
% end
% K_array = [8.063967e+04; 6.026747e+05 ; 5.531480e+05 ; 6.425025e+05;6.258761e+05;  1.770501e+05];
K_initial=build_stiffness_matrix(k_opt);
[Phi,D] = eig(K_initial,M);

lambda = diag(D);

% Remove numerical errors
lambda(lambda<0)=0;

omega = sqrt(lambda);

frequency = omega/(2*pi);

%% ---------------------------------------------------------
% SORT MODES
% ----------------------------------------------------------

[frequency,index] = sort(frequency);

Phi = Phi(:,index);

%% ---------------------------------------------------------
% PRINT NATURAL FREQUENCIES
% ----------------------------------------------------------

disp(' ');
disp('Natural Frequencies (Hz)');
disp('-------------------------');

for i=1:n

    fprintf('Mode %d : %.4f Hz\n',i,frequency(i));

end

%% ---------------------------------------------------------
% NORMALIZE MODE SHAPES
% ----------------------------------------------------------

for i=1:n

    Phi(:,i)=Phi(:,i)/max(abs(Phi(:,i)));

end

%% ---------------------------------------------------------
% PLOT MODE SHAPES
% ----------------------------------------------------------

stories = 1:n;

figure('Color','w');

hold on

colors = lines(n);

for i=1:n

    plot(Phi(:,i),stories,...
        '-o',...
        'LineWidth',2,...
        'MarkerSize',8,...
        'Color',colors(i,:),...
        'DisplayName',['Mode ',num2str(i)]);

end

xlabel('Normalized Mode Shape','FontSize',12)
ylabel('Storey','FontSize',12)
title('6-DOF Shear Building Mode Shapes','FontSize',14)

grid on
box on

legend('Location','best')

set(gca,...
    'FontSize',12,...
    'LineWidth',1.2,...
    'YDir','normal')
