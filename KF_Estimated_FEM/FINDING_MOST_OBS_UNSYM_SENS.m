a = 0.6; b = 0.4;
num_sensors = 8;

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
title('8 Optimized Asymmetric Sensors for 6 Modes');
axis([0 a 0 b]);
%% SPECIFICALLY THIS CONFIGURATION 
x_sens = [0.1292, 0.1402 , 0.1402 , 0.2945 , 0.3055, 0.3055 ,0.4598 ,0.4708];
y_sens = [0.2037, 0.0788,0.3212,0.1963,  0.0861,0.3139,0.0788,  0.1963];

% Mode Indices
m = [1, 1, 2, 1, 3, 3];
n = [1, 2, 1, 3, 1, 3];

Phi = zeros(8,6);
for i = 1:8
    for j = 1:6
        Phi(i,j) = sin(m(j)*pi*x_sens(i)/a) * sin(n(j)*pi*y_sens(i)/b);
    end
end

fprintf('Rank: %d\n', rank(Phi)); 
fprintf('Condition Number: %.4f\n', cond(Phi));
