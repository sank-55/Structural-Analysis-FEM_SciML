% Free Vibration of Beam using FEM
clear; clc; close all

% Beam properties
E = 210e9;         % Young's modulus (Pa)
rho = 8050;        % Density (kg/m^3)
A = 40e-6;         % Cross-sectional area (m^2)
I = 1e-6;          % Moment of inertia (m^4)
L_total = 1;       % Total length (m)

% FEM Discretization
nElem = 40;                  % Number of elements
nNode = nElem + 1;            % Number of nodes
dx = L_total / nElem;         % Length of each element

% Global matrices initialization
K_global = zeros(2*nNode);    % Stiffness matrix
M_global = zeros(2*nNode);    % Mass matrix

% Element matrices
for e = 1:nElem
    % Local stiffness matrix (Euler-Bernoulli beam)
    Ke = E*I/dx^3 * [12,    6*dx,   -12,    6*dx;
                     6*dx,  4*dx^2, -6*dx,  2*dx^2;
                    -12,   -6*dx,    12,   -6*dx;
                     6*dx,  2*dx^2, -6*dx,  4*dx^2];
    
    % Local mass matrix (consistent mass matrix)
    m = rho*A*dx;
    Me = m/420 * [156,   22*dx,   54,    -13*dx;
                  22*dx, 4*dx^2, 13*dx,  -3*dx^2;
                  54,    13*dx,  156,    -22*dx;
                 -13*dx, -3*dx^2, -22*dx, 4*dx^2];

    % Assembly
    idx = [2*e-1, 2*e, 2*e+1, 2*e+2];
    K_global(idx, idx) = K_global(idx, idx) + Ke;
    M_global(idx, idx) = M_global(idx, idx) + Me;
end

% Apply boundary conditions (Simply Supported Beam at both ends)
fixed_dofs = [1, 2, 2*nNode-1, 2*nNode]; % Remove displacement and slope at both ends

free_dofs = setdiff(1:2*nNode, fixed_dofs);

K_reduced = K_global(free_dofs, free_dofs);
M_reduced = M_global(free_dofs, free_dofs);

% Solve eigenvalue problem
[phi, omega2] = eig(K_reduced, M_reduced);

% Natural frequencies (Hz)
omega = sqrt(diag(omega2));       % rad/s
freq = omega / (2*pi);             % Hz

% Sort frequencies
[freq_sorted, idx_sort] = sort(freq);
phi = phi(:, idx_sort);

% Display first few natural frequencies
disp('First few natural frequencies (Hz):');
disp(freq_sorted(1:5));

% Plot first 3 mode shapes
nModes = 3;
x = linspace(0, L_total, nNode);

for i = 1:nModes
    mode_shape = zeros(2*nNode,1);
    mode_shape(free_dofs) = phi(:,i);

    figure;
    plot(x, mode_shape(1:2:end), '-o', 'LineWidth', 2);
    title(['Mode Shape ', num2str(i), ' - f = ', num2str(freq_sorted(i), '%.2f'), ' Hz']);
    xlabel('Beam Length (m)');
    ylabel('Transverse Displacement (arbitrary)');
    grid on;
end
