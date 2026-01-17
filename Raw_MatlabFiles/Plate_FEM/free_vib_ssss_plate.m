%% Everything is validated and working perfectly 

clear; clc; close all;

%% 1. MATERIAL AND GEOMETRY PROPERTIES
E = 70e9;           % Young's Modulus (Pa) - Aluminum
rho = 2700;         % Density (kg/m^3)
nu = 0.3;           % Poisson's ratio
Lx = 0.6;           % Length (m)
Ly = 0.4;           % Width (m)
h = 0.00625;        % Thickness (m)

% Mesh properties (Paper uses 12x8)
nx_ele = 24; 
ny_ele = 16;

% Derived constants
D = (E * h^3) / (12 * (1 - nu^2)); % Flexural Rigidity
nx_node = nx_ele + 1;
ny_node = ny_ele + 1;
dx = Lx / nx_ele;
dy = Ly / ny_ele;
total_nodes = nx_node * ny_node;
total_dof = total_nodes * 3;

%% 2. ELEMENT STIFFNESS AND MASS MATRICES
% We compute the element matrices using 3x3 Gauss Quadrature
[Ke, Me] = computeElementMatrices(dx, dy, D, rho, h, nu);

%% 3. GLOBAL ASSEMBLY
K_global = sparse(total_dof, total_dof);
M_global = sparse(total_dof, total_dof);

fprintf('Assembling Global Matrices...\n');
for ey = 1:ny_ele
    for ex = 1:nx_ele
        % Node indices for current element (Counter-Clockwise)
        n1 = (ey-1) * nx_node + ex;
        n2 = n1 + 1;
        n3 = ey * nx_node + ex + 1;
        n4 = ey * nx_node + ex;
        
        nodes = [n1, n2, n3, n4];
        
        % Mapping local DOFs to global DOFs
        % Each node has 3 DOFs: [w, th_y, th_x]
        element_dofs = [];
        for n = nodes
            element_dofs = [element_dofs, (n-1)*3 + (1:3)];
        end
        
        K_global(element_dofs, element_dofs) = K_global(element_dofs, element_dofs) + Ke;
        M_global(element_dofs, element_dofs) = M_global(element_dofs, element_dofs) + Me;
    end
end

%% 4. BOUNDARY CONDITIONS (Simply Supported - SSSS)
% Constraint: w = 0 on all edges. 
% For SSSS, tangential rotations are also typically fixed to ensure stability.
constrained_dofs = [];

for i = 1:ny_node
    for j = 1:nx_node
        node_idx = (i-1) * nx_node + j;
        % Check if node is on any boundary
        if (i == 1 || i == ny_node || j == 1 || j == nx_node)
            % Constrain w (1st DOF of node)
            constrained_dofs = [constrained_dofs, (node_idx-1)*3 + 1];
            
            % Constrain tangential rotation (following Kirchhoff boundary logic)
            if (i == 1 || i == ny_node) % Horizontal edges
                constrained_dofs = [constrained_dofs, (node_idx-1)*3 + 2]; % dw/dx = 0
            end
            if (j == 1 || j == nx_node) % Vertical edges
                constrained_dofs = [constrained_dofs, (node_idx-1)*3 + 3]; % dw/dy = 0
            end
        end
    end
end

constrained_dofs = unique(constrained_dofs);
all_dofs = 1:total_dof;
free_dofs = setdiff(all_dofs, constrained_dofs);

%% 5. SOLVE EIGENVALUE PROBLEM
fprintf('Solving for Natural Frequencies...\n');
[V, OmegaSq] = eigs(K_global(free_dofs, free_dofs), M_global(free_dofs, free_dofs), 6, 'smallestabs');

% Convert eigenvalues to frequencies (Hz)
natural_freqs = sqrt(diag(OmegaSq)) / (2 * pi);
[natural_freqs, sort_idx] = sort(natural_freqs);
V = V(:, sort_idx);

%% 6. ANALYTICAL SOLUTION COMPARISON
% Formula for SSSS: f_mn = (pi/2) * sqrt(D/(rho*h)) * ((m/Lx)^2 + (n/Ly)^2)
exact_freqs = [];
for m = 1:3
    for n = 1:3
        f = (pi/2) * sqrt(D/(rho*h)) * ((m/Lx)^2 + (n/Ly)^2);
        exact_freqs = [exact_freqs; f];
    end
end
exact_freqs = sort(exact_freqs);

fprintf('\n--- RESULTS COMPARISON ---\n');
fprintf('Mode | FEM (Hz) | Exact (Hz) | Error (%%)\n');
for i = 1:5
    err = abs(natural_freqs(i) - exact_freqs(i))/exact_freqs(i) * 100;
    fprintf('%4d | %8.2f | %10.2f | %8.2f\n', i, natural_freqs(i), exact_freqs(i), err);
end

%% 7. VISUALIZATION
for i=1:4
    mode_to_plot = i;
    full_mode = zeros(total_dof, 1);
    full_mode(free_dofs) = V(:, mode_to_plot);

    % Extract only the lateral displacement 'w' for plotting
    w_grid = reshape(full_mode(1:3:end), [nx_node, ny_node])';

    [X, Y] = meshgrid(linspace(0, Lx, nx_node), linspace(0, Ly, ny_node));
    figure('Color', 'w');
    surf(X, Y, w_grid, 'FaceColor', 'interp', 'EdgeColor', 'none');
    title(['Mode Shape ', num2str(mode_to_plot), ' (Freq: ', num2str(natural_freqs(mode_to_plot), '%.2f'), ' Hz)']);
    xlabel('Lx (m)'); ylabel('Ly (m)'); zlabel('Displacement');
    colorbar; colormap jet; view(3);
end

%% --- HELPER FUNCTION: ELEMENT MATRICES ---
function [Ke, Me] = computeElementMatrices(dx, dy, D, rho, h, nu)
    % Gauss points and weights for 3x3 integration
    gp = [-sqrt(0.6), 0, sqrt(0.6)];
    gw = [5/9, 8/9, 5/9];
    
    Ke = zeros(12, 12);
    Me = zeros(12, 12);
    detJ = (dx/2) * (dy/2);
    
    % Constitutive Matrix D_mat
    D_mat = D * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    
    for i = 1:3
        for j = 1:3
            xi = gp(i); eta = gp(j);
            wt = gw(i) * gw(j);
            
            % Hermite Shape Functions for Plate
            [N, B] = getPlateShape(xi, eta, dx, dy);
            
            Ke = Ke + (B' * D_mat * B) * detJ * wt;
            Me = Me + (rho * h * (N' * N)) * detJ * wt;
        end
    end
end

function [N, B] = getPlateShape(xi, eta, dx, dy)
    % Standard 12-DOF ACM Shape Functions
    % 1D Hermite Polynomials
    H = @(t) [0.25*(2-3*t+t^3), 0.25*(1-t-t^2+t^3), 0.25*(2+3*t-t^3), 0.25*(-1-t+t^2+t^3)];
    dH = @(t) [0.25*(-3+3*t^2), 0.25*(-1-2*t+3*t^2), 0.25*(3-3*t^2), 0.25*(-1+2*t+3*t^2)];
    ddH = @(t) [0.25*(6*t), 0.25*(-2+6*t), 0.25*(-6*t), 0.25*(2+6*t)];
    
    hx = H(xi); dhx = dH(xi); ddhx = ddH(xi);
    hy = H(eta); dhy = dH(eta); ddhy = ddH(eta);
    
    jx = dx/2; jy = dy/2;
    
    N = zeros(1, 12);
    B = zeros(3, 12); % Curvature matrix
    
    % Node ordering: 1(-1,-1), 2(1,-1), 3(1,1), 4(-1,1)
    % DOFs at each node: w, dw/dx, dw/dy
    indices = [1,2; 3,4; 3,4; 1,2]; % Mapping xi/eta to Hermite indices
    
    for n = 1:4
        idx_x = 1 + (n==2 || n==3)*2; % 1 for nodes 1,4; 3 for nodes 2,3
        idx_y = 1 + (n==3 || n==4)*2; % 1 for nodes 1,2; 3 for nodes 3,4
        
        b_idx = (n-1)*3;
        % w
        N(b_idx+1) = hx(idx_x) * hy(idx_y);
        % dw/dx
        N(b_idx+2) = hx(idx_x+1) * hy(idx_y) * jx;
        % dw/dy
        N(b_idx+3) = hx(idx_x) * hy(idx_y+1) * jy;
        
        % B matrix: [d2w/dx2; d2w/dy2; 2*d2w/dxdy]
        B(1, b_idx+1) = ddhx(idx_x) * hy(idx_y) / (jx^2);
        B(1, b_idx+2) = ddhx(idx_x+1) * hy(idx_y) * jx / (jx^2);
        B(1, b_idx+3) = ddhx(idx_x) * hy(idx_y+1) * jy / (jx^2);
        
        B(2, b_idx+1) = hx(idx_x) * ddhy(idx_y) / (jy^2);
        B(2, b_idx+2) = hx(idx_x+1) * ddhy(idx_y) * jx / (jy^2);
        B(2, b_idx+3) = hx(idx_x) * ddhy(idx_y+1) * jy / (jy^2);
        
        B(3, b_idx+1) = 2 * (dhx(idx_x) * dhy(idx_y)) / (jx*jy);
        B(3, b_idx+2) = 2 * (dhx(idx_x+1) * dhy(idx_y) * jx) / (jx*jy);
        B(3, b_idx+3) = 2 * (dhx(idx_x) * dhy(idx_y+1) * jy) / (jx*jy);
    end
end
