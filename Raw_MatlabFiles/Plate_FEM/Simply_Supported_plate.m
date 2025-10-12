clear; close all;
%%%%% __________ BENDING OF THIN_PLATES PROBLEM with more Refinement of the mesh __________________________%%%%
% where w = w1N1 + w2N2 +....+ w12N12
%isoparametric form
%% Properties of the plate 
% Inputs:
E = 210e9;                      % Young's modulus (Pa)
h = 0.005;                     %  Thickness of the plate in m
rho = 8050;                     % density of the plate material(kg/m3)
q= -1500;                           % Transvere (perpendicular) distributed force over domain (N/m2)
nu = 0.3;                       % poisons ratio
le = 1;                         % elemental length 
D = (E*h^3 / (12*(1 - nu^2))) * [1, nu, 0;
                                 nu, 1, 0;
                                 0, 0, (1 - nu)/2];
% put Transverse applied load( in N)
% right edge ; % left edge
%Vx =5000 ;      Vy =12500;

% Plate dimensions
Lx = 2; % Length along x
Ly = 2; % Length along y

% Number of elements along x and y ( Define such a way so that dx=dy)
nx = 20;   % along x
ny = 20;   % along y
dx = Lx/nx; % elemntal length in x
dy = Ly/ny; % elemntal length in y
Nelem = nx * ny;       % Total number of elements
Nnode = (nx+1)*(ny+1); % Total number of nodes
Ndof = 3*Nnode;        % 3 dofs per node

%% Generate nodal coordinates
x = linspace(0, Lx, nx+1);
y = linspace(0, Ly, ny+1);
%[X, Y] = meshgrid(x, y);
% Generate node coordinates
coords = zeros(Nnode, 2);
node = 1;
for j = 0:ny
    for i = 0:nx
        coords(node, :) = [i*dx, j*dy];
        node = node + 1;
    end
end

%% Generate connectivity                  % in a CCW manner 
c = zeros(Nelem, 4);                        %  .  .  .  . ...
elem = 0;                                   %  5--6--7--8....
% logic for connectivity matrix             %  |  |  |  |
for j = 1:ny                                %  1--2--3--4....
    for i = 1:nx
        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n3 = n2 + (nx+1);
        n4 = n1 + (nx+1);
        elem = elem + 1;
        c(elem, :) = [n1, n2, n3, n4];
    end
end

%% Gauss points(2D) and weights for 2x2 integration
gp = [-1/sqrt(3) ,  1/sqrt(3)];
w= [1 , 1];

%% Initialize the matrix
Kg = zeros(Ndof, Ndof);
Mg = zeros(Ndof, Ndof);
Fg = zeros(Ndof,1); 

%%%%%% ********    % SHAPE FUNCTION & ITS DERIVATIVES %%%%%%%%%%
syms xi eta real 
           % 1D Hermite shape functions (cubic), defined over [-1, 1]
        H1 = 1 - 3*xi^2 + 2*xi^3;           % Displacement at node i
        H2 = le*(xi - 2*xi^2 + xi^3);        % Slope at node i
        H3 = 3*xi^2 - 2*xi^3;               % Displacement at node j
         H4 = le*(-xi^2 + xi^3);              % Slope at node j

        H_xi  = [H1; H2; H3; H4];   % Functions in xi direction
        %H = [H1; H4; H3; H2];
        H_eta = subs(H_xi, xi, eta);        % Functions in eta direction

        %% Shape Function (Neglecting Shear)
            N1  = H_xi(1) * H_eta(1);  % w at Node 1
            N2  = H_xi(2) * H_eta(1);  % θx at Node 1
            N3  = H_xi(1) * H_eta(2);  % θy at Node 1
            N4  = H_xi(3) * H_eta(1);  % w at Node 2
            N5  = H_xi(4) * H_eta(1);  % θx at Node 2
            N6  = H_xi(3) * H_eta(2);  % θy at Node 2
            N7  = H_xi(3) * H_eta(3);  % w at Node 3
            N8  = H_xi(4) * H_eta(3);  % θx at Node 3
            N9  = H_xi(3) * H_eta(4);  % θy at Node 3
            N10  = H_xi(1) * H_eta(3);  % w at Node 4
            N11  = H_xi(2) * H_eta(3);  % θx at Node 4
            N12  = H_xi(1) * H_eta(4);  % θy at Node 4

    %SHAPE FUNCTION IN MATRIX FORM
    N=[N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12];
    
    %% Compute 1st Derivatives of N w.r.t. xi and eta
    N_xi  = simplify(diff(N, xi));       % ∂N/∂xi
    N_eta = simplify(diff(N, eta));     % ∂N/∂eta
    dN_dxi =[ N_xi;
              N_eta];
    %% Compute 2nd Derivatives of N w.r.t. xi and eta

    N_xixi  = simplify(diff(N, xi, 2));       % ∂²N/∂xi²
    N_etaeta = simplify(diff(N, eta, 2));     % ∂²N/∂eta²
    N_xieta = simplify(diff(diff(N, xi), eta)); % ∂²N/∂xi∂eta

%% Computation (HAVE TO LOOK IN THIS PART ))
for e = 1:Nelem
    nodes = c(e,:);
    XY = coords(nodes,:); % takes the 4 coordinares of perticular element
    % Derivatives w.r.t. xi and eta
    dN_dxi = zeros(2,12);
    %local stiffness matrix initialization within a element
    Ke = zeros(12,12);
    Me = zeros(12,12);
    Fe = zeros(12,1);
    for i = 1:2
        for j = 1:2
            xi_val = gp(i);
            eta_val = gp(j);
            wts = w(i) * w(j);

            % Symbolic substitution of xi and eta
            N_val = double(subs(N, [xi, eta], [xi_val, eta_val]));
            N_xixi_val = double(subs(N_xixi, [xi, eta], [xi_val, eta_val]));
            N_etaeta_val = double(subs(N_etaeta, [xi, eta], [xi_val, eta_val]));
            N_xieta_val = double(subs(N_xieta, [xi, eta], [xi_val, eta_val]));
           %%%%%% *************Jacobian it may change per one mapping and element wise
           % FOR linear map xi--x & eta--y (in this case dx=dy)
           J = [2/dx 0;
                0  2/dy];
           %% J = dN_dxi* XY;     %           
            detJ = det(J);                     % Jacobian determinant
            %dN_dx = inv(J)*dN_dxi;              %basically dN_dx = [dN/dx, dN/dy] matrix

            % ---- Calculating B MAtrix  ----
            B = zeros(3, 12);
            % %*** Will use in case of nonlinear mapping
            % V =zeros(5, 12);
            %  % assignning the B matrix change according to element
            %        V = [ N_xixi;
            %            N_etaeta;
            %            N_xieta;
            %            N_xi;
            %            N_eta];
            % 

                 B =  (4/dx^2)*[ N_xixi_val;    % [d2N1/dx2   d2N2/dx2 ..... d2Nn/dx2 
                                 N_etaeta_val;  %  d2N1/dy2   d2N2/dy2 ..... d2Nn/dy2 
                                2*N_xieta_val];  %  d2N1/dxdy  d2N2/dxdy ..... d2Nn/dxdy]

            % Assigning the values into elemental local matrix
            Ke = Ke + B' * D * B * detJ * wts;  %% though For linear mapping le not going to affect the ke 
            Me = Me + ((rho*h)*(detJ*wts))*(N_val')*N_val;
            % Calculating global force vector
            Fe = Fe +  q*detJ*wts*(N_val')  ;
        end
       %Fe = Fe +  q*detJ*wts*N_val'  ;
    end

    %%%% ***Put into the logic for finding global k
    dofs = reshape([3*nodes-2; 3*nodes-1; 3*nodes],1,[]);
    % Global matrix is properly shaped before adding ke in required position
    % Add to global stiffness matrix
    Kg(dofs,dofs) = Kg(dofs,dofs) + Ke; %% **reshape funcn
    Mg(dofs,dofs) = Mg(dofs,dofs) + Me;
    % Global force vector computing ****
    Fg(dofs,1) = Fg(dofs,1) + Fe;
end


%% Boundary conditions:
% ****put varrious combination****
%  nodal dof ( w , dw/dx and dw/dy ) for each of the element

% % Find nodes on left edge (x=0)
% left_nodes = find(abs(coords(:,1)) < 1e-6);
% 
% % Their DOFs
% fixedDOF = reshape([3*left_nodes-2; 3*left_nodes-1; 3*left_nodes],1,[]);                           %---(% Write what is the boundary condn)--
% 
% 
% % set the rest of the dof as free
% freeDOF = setdiff(1:Ndof, fixedDOF); % setdiff remove all the fixed nodes from the dof vector

% Find nodes on all 4 edges
tol = 1e-8;   % tolerance for floating point comparison
left_nodes   = find(abs(coords(:,1)) < tol);
right_nodes  = find(abs(coords(:,1) - Lx) < tol);
bottom_nodes = find(abs(coords(:,2)) < tol);
top_nodes    = find(abs(coords(:,2) - Ly) < tol);

% All edge nodes (no duplicates)
edge_nodes = unique([left_nodes; right_nodes; bottom_nodes; top_nodes]);

% Clamped DOFs: for each node, w, theta_x, theta_y must be zero
% fixedDOF = reshape([3*left_nodes-2; 3*left_nodes-1; 3*left_nodes],1,[]);
% 
% Free DOFs
% freeDOF = setdiff(1:Ndof, fixedDOF);

% Simply supported DOFs: only w (transverse displacement) is zero on edges
fixedDOF = 3 * edge_nodes - 2;   % w DOFs only

% Free DOFs
freeDOF = setdiff(1:Ndof, fixedDOF);




% Apply boundary 
Kf = Kg(freeDOF, freeDOF); %%fianl stiffness matrix
Mf = Mg(freeDOF, freeDOF); %%fianl Mass-coeff matrix
%   assign the load after understanding ## use the reshape to get the
%   reshape the Fg then assign the external load
Ff = Fg(freeDOF) ;%+[0;0;0;0;0;0] ;

%% ---- Solving the system ----
%calculation of the deflection and slopes
U = zeros(Ndof, 1);
U(freeDOF) = Kf \ Ff;

% for finding the force
F_f= Kg*U  ; %final force vector 


%% Calculating eig values  ------
    % Solve eigenvalue problem
[phi, omega2] = eig(Kf, Mf);

    % Natural frequencies (Hz)
omega = sqrt(diag(omega2));       % rad/s
freq = omega / (2*pi);             % Hz

    % Sorting frequencies
[freq_sorted, idx_sort] = sort(freq);
phi = phi(:, idx_sort);

% Display first few natural frequencies
disp('First few natural frequencies (Hz):');
disp(freq_sorted(1:5));



%% Plotting the Modeshapes of flat plate
nModes = 4;  % Number of modes to plot

% Sort nodes to improve plotting grid
[~, idx] = sortrows(coords, [2 1]);
coordsS = coords(idx, :);  % sorted coordinates

% High-resolution grid for interpolation
xRange = linspace(min(coords(:,1)), max(coords(:,1)), 150);
yRange = linspace(min(coords(:,2)), max(coords(:,2)), 150);
[Xq, Yq] = meshgrid(xRange, yRange);  % query grid

for m = 1:nModes
    % Construct full mode shape
    modeShape = zeros(Ndof, 1);
    modeShape(freeDOF) = phi(:, m);
    w = modeShape(1:3:end);         % Only transverse displacement

    % Interpolate to smooth surface
    Wq = griddata(coords(:,1), coords(:,2), w, Xq, Yq, 'cubic');

    % Plotting
    figure;
    surf(Xq, Yq, Wq, 'EdgeColor', 'none');  % no mesh lines
    shading interp;
    colormap turbo;
    colorbar;
    axis equal tight;
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('w (m)');
    title(sprintf('Mode %d, f = %.2f Hz', m, freq_sorted(m)));
end


%%   ==== 3D DEFORMATION PLOT OF PLATE ====
% Extract real physical transverse deflection (no normalization)
W = U(1:3:end);   % w DOF from [w, theta_x, theta_y]

% Coordinates
X = coords(:,1);
Y = coords(:,2);

% Triangulate mesh (for unstructured or Q4 grid)
tri = delaunay(X, Y);

% OPTIONAL: Visual scale factor for visibility (can set to 1 for exact)
scaleFactor = 1.0;  % Set to e.g. 10 to exaggerate, or 1 for real scale
Wvis = scaleFactor * W;

% --- Plotting physically-realistic 3D deformation ---
figure('Color','w');

trisurf(tri, X, Y, Wvis, ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'SpecularStrength', 0.3);

shading interp;
colormap turbo;
colorbar;

xlabel('X [m]', 'FontSize', 12);
ylabel('Y [m]', 'FontSize', 12);
zlabel('Deflection w [m]', 'FontSize', 12);
title(sprintf('Realistic Plate Deflection (Kirchhoff FEM)\nMax w = %.4f m', max(abs(W))), ...
      'FontSize', 14, 'FontWeight','bold');

% 3D view and lighting
view([-35 25]); axis equal tight;
camlight('headlight'); lighting gouraud;
material shiny;


%% 2d subplot
% ==== Extract DOFs ====
w_vals  = U(1:3:end);    % vertical deflection
tx_vals = U(2:3:end);    % slope along x (theta_x)
ty_vals = U(3:3:end);    % slope along y (theta_y)

X = coords(:,1);
Y = coords(:,2);

% ==== Midline in X-direction (y = mid-height) ====
y_mid = mean(unique(Y));  % Midline in Y-direction (y = constant)
tol = 1e-5;                % Tolerance to select midline nodes
x_line_idx = find(abs(Y - y_mid) < tol);  % Select nodes along midline
[X_xline, sortX] = sort(X(x_line_idx));   % Sort X positions
w_xline = w_vals(x_line_idx(sortX));     % Corresponding deflections along X
tx_xline = tx_vals(x_line_idx(sortX));   % Corresponding slopes along X
ty_xline = ty_vals(x_line_idx(sortX));   % Corresponding slopes along X

% Smooth tx and ty using spline interpolation for smoothness
tx_xline_smooth = spline(X_xline, tx_xline, X_xline); 
ty_xline_smooth = spline(X_xline, ty_xline, X_xline);

% ==== Midline in Y-direction (x = mid-width) ====
x_mid = mean(unique(X));  % Midline in X-direction (x = constant)
y_line_idx = find(abs(X - x_mid) < tol);  % Select nodes along midline
[Y_yline, sortY] = sort(Y(y_line_idx));   % Sort Y positions
w_yline = w_vals(y_line_idx(sortY));     % Corresponding deflections along Y
tx_yline = tx_vals(y_line_idx(sortY));   % Corresponding slopes along Y
ty_yline = ty_vals(y_line_idx(sortY));   % Corresponding slopes along Y

% Smooth tx and ty using spline interpolation for smoothness
tx_yline_smooth = spline(Y_yline, tx_yline, Y_yline);
ty_yline_smooth = spline(Y_yline, ty_yline, Y_yline);

% ==== Plot Variation Along X ====
figure;
subplot(3,1,1)
plot(X_xline, w_xline, '-o','LineWidth',1.5);
ylabel('Deflection w [m]'); 
title('Variation Along x-direction (y = mid)');
grid on;

subplot(3,1,2)
plot(X_xline, tx_xline_smooth, '-o','LineWidth',1.5);
ylabel('\theta_x [rad]');
title('Variation of \theta_x along x');
grid on;

subplot(3,1,3)
plot(X_xline, ty_xline_smooth, '-o','LineWidth',1.5);
ylabel('\theta_y [rad]'); 
xlabel('x [m]');
title('Variation of \theta_y along x');
grid on;

% ==== Plot Variation Along Y ====
figure;
subplot(3,1,1)
plot(Y_yline, w_yline, '-o','LineWidth',1.5);
ylabel('Deflection w [m]'); 
title('Variation Along y-direction (x = mid)');
grid on;

subplot(3,1,2)
plot(Y_yline, tx_yline_smooth, '-o','LineWidth',1.5);
ylabel('\theta_x [rad]');
title('Variation of \theta_x along y');
grid on;

subplot(3,1,3)
plot(Y_yline, ty_yline_smooth, '-o','LineWidth',1.5);
ylabel('\theta_y [rad]'); 
xlabel('y [m]');
title('Variation of \theta_y along y');
grid on;



%% Plotting undeforemed mesh
figure;
hold on; axis equal;
title('Undeformed Mesh');
xlabel('X'); ylabel('Y');

% Loop through each element and draw its edges
for e = 1:size(c,1)
    nodes = c(e,:);        % Node indices of the element
    xe = coords(nodes,1);      % X coordinates of element nodes
    ye = coords(nodes,2);      % Y coordinates of element nodes
    
    % Close the loop by appending the first node at the end
    patch(xe([1:end 1]), ye([1:end 1]), 'w', 'EdgeColor', 'k');
end

hold off;
