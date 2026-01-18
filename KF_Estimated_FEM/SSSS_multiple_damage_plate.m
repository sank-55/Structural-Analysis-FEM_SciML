%% Plate using dual kalman filter 


%% ---------------------------  There is lot more to perfectly tune the P,Q,R quantities and understand the undergoing physics ----


clear; close all; clc;
syms x y

%% Given parameters 
a=0.6;
b=0.4;
q = -100; % N/m2 

xdiv=12;
ydiv=8;

Lex=a/xdiv;
Ley=b/ydiv;
rho=2700;
E=70e9;
h=0.00625;
neu=0.3;
D = E*h^3/(12*(1-neu^2));



dof_per_node=4; % dof per node 
Ne=(xdiv)*(ydiv); % total no. of element 
Nnode=(xdiv+1)*(ydiv+1); % total no. of the nodes 
Ndof = Nnode*dof_per_node;

%% Calculating the coordinates 
coords = zeros(Nnode,2);
n = 1;
for j = 1:ydiv+1
    for i = 1:xdiv+1
        coords(n, :) = [(i-1)*Lex, (j-1)*Ley];
        % node_num(j, i) = n;
        n = n + 1;
    end
end



%% applying Boundary condition 
BCdof=[1:4:(xdiv+1)*dof_per_node     ydiv*(xdiv+1)*dof_per_node+1:4:(ydiv+1)*(xdiv+1)*dof_per_node    2:4:(xdiv+1)*dof_per_node       ydiv*(xdiv+1)*dof_per_node+2:4:(ydiv+1)*(xdiv+1)*dof_per_node];


count=4*(xdiv+1);
for i=1:ydiv-1
    count=count+1;
    BCdof(count)=(xdiv+1)*dof_per_node+(i-1)*(xdiv+1)*dof_per_node+1;
%     count=count+1;
%     BCdof(count)=(xdiv+1)*Nedof+(i-1)*(xdiv+1)*Nedof+2;
    count=count+1;
    BCdof(count)=(xdiv+1)*dof_per_node+(i-1)*(xdiv+1)*dof_per_node+3;
% %     count=count+1;
%     BCdof(count)=(xdiv+1)*Nedof+(i-1)*(xdiv+1)*Nedof+4;

    count=count+1;
    BCdof(count)=(xdiv+1)*dof_per_node+(i)*(xdiv+1)*dof_per_node-3;
%     count=count+1;
%     BCdof(count)=(xdiv+1)*Nedof+(i)*(xdiv+1)*Nedof-2;
    count=count+1;
     BCdof(count)=(xdiv+1)*dof_per_node+(i)*(xdiv+1)*dof_per_node-1;
% %     count=count+1;
%     BCdof(count)=(xdiv+1)*Nedof+(i)*(xdiv+1)*Nedof;
end

Fdof=setdiff([1:Ndof],BCdof);


%% Shape or Interpolation functions 
Nx=[(2*x^3-3*x^2*Lex+Lex^3)/Lex^3 (x^3*Lex-2*x^2*Lex^2+x*Lex^3)/Lex^3 (-2*x^3+3*x^2*Lex)/Lex^3 (x^3*Lex-x^2*Lex^2)/Lex^3];
Ny=[(2*y^3-3*y^2*Ley+Ley^3)/Ley^3 (y^3*Ley-2*y^2*Ley^2+y*Ley^3)/Ley^3 (-2*y^3+3*y^2*Ley)/Ley^3 (y^3*Ley-y^2*Ley^2)/Ley^3];

% creating the N matrix 
N(1:4)=[Nx(1)*Ny(1) Nx(2)*Ny(1) Nx(1)*Ny(2) Nx(2)*Ny(2)];
N(5:8)=[Nx(3)*Ny(1) Nx(4)*Ny(1) Nx(3)*Ny(2) Nx(4)*Ny(2)];
N(9:12)=[Nx(3)*Ny(3) Nx(4)*Ny(3) Nx(3)*Ny(4) Nx(4)*Ny(4)];
N(13:16)=[Nx(1)*Ny(3) Nx(2)*Ny(3) Nx(1)*Ny(4) Nx(2)*Ny(4)];
% taking the differentiation with respect to x,y 
Nxx=diff(N,x,2);
Nyy=diff(N,y,2);
Nxy=diff(diff(N,x),y);

%% The elemental mass and stiffness function using symbollic intergration
Me=rho*h*double(int(int(N'*N,x,0,Lex),y,0,Ley));

Ke=D*double(int(int((Nxx'*Nxx+2*Nxy'*Nxy+Nyy'*Nyy),x,0,Lex),y,0,Ley));

Fe = q*double(int(int(N',x,0,Lex),y,0,Ley));
%% Global assembly 
% FOR DAMAGED MATRIX
Ke_store = cell(Ne,1); % store element stiffness
elemDof = cell(Ne,1);
conn_nodes = zeros(Ne,4); % as we using 4 noded element gives Connectivity of the nodes
elemCenter = zeros(Ne,2 );

count=0;

for i=1:xdiv
    for j=1:ydiv
        count=count+1;
        C(count,1:4)=[dof_per_node*(j-1)*(xdiv+1)+(i-1)*dof_per_node+1 dof_per_node*(j-1)*(xdiv+1)+(i-1)*dof_per_node+2 dof_per_node*(j-1)*(xdiv+1)+(i-1)*dof_per_node+3 dof_per_node*(j-1)*(xdiv+1)+(i-1)*dof_per_node+4];
        conn_nodes(count,1) = 1 + (C(count,1)-1)/4 ;
        C(count,5:8)=[dof_per_node*(j-1)*(xdiv+1)+(i-1)*dof_per_node+5 dof_per_node*(j-1)*(xdiv+1)+(i-1)*dof_per_node+6 dof_per_node*(j-1)*(xdiv+1)+(i-1)*dof_per_node+7 dof_per_node*(j-1)*(xdiv+1)+(i-1)*dof_per_node+8];
        conn_nodes(count,2) = 1 + (C(count,5)-1)/4 ;
        C(count,13:16)=[dof_per_node*(j)*(xdiv+1)+(i-1)*dof_per_node+1 dof_per_node*(j)*(xdiv+1)+(i-1)*dof_per_node+2 dof_per_node*(j)*(xdiv+1)+(i-1)*dof_per_node+3 dof_per_node*(j)*(xdiv+1)+(i-1)*dof_per_node+4];
        conn_nodes(count,4) = 1 + (C(count,13)-1)/4 ;
        C(count,9:12)=[dof_per_node*(j)*(xdiv+1)+(i-1)*dof_per_node+5 dof_per_node*(j)*(xdiv+1)+(i-1)*dof_per_node+6 dof_per_node*(j)*(xdiv+1)+(i-1)*dof_per_node+7 dof_per_node*(j)*(xdiv+1)+(i-1)*dof_per_node+8];
        conn_nodes(count,3) = 1 + (C(count,9)-1)/4 ;
    end
end

K=zeros(dof_per_node*Nnode,dof_per_node*Nnode);
M=zeros(dof_per_node*Nnode,dof_per_node*Nnode);
F = zeros(Ndof,1);
for e=1:Ne
    nodes = conn_nodes(e, :);
    elemCenter(e,:) = mean(coords(nodes,:),1); % store the centroids 
     for i=1:dof_per_node*4
         for j=1:dof_per_node*4
             M(C(e,i),C(e,j))=M(C(e,i),C(e,j))+Me(i,j);
             K(C(e,i),C(e,j))=K(C(e,i),C(e,j))+Ke(i,j);
             
         end
         F(C(e,i),1) = F(C(e,i),1) + Fe(i,1);
     end
     elemDof{e} = C(e,:);
     Ke_store{e} = Ke; % for damage matrix 
end
       


%% =======  ZONING THE DOMAIN ==========
% FOR This case only one zone is taken
Nz = 6 ;
zoneID = sparse(Ne,1);


%% === divided into several zones ====================
for e = 1:Ne
    x = elemCenter(e,1);
    y= elemCenter(e,2);
    if x <= a/3 
        if y <= b/2
            zoneID(e) = 1;
        else 
            zoneID(e) = 4;
        end
    elseif x <= 2*a/3
        if y <= b/2 
            zoneID(e) = 2;
        else
            zoneID(e) = 5;
        end
    elseif x <= 3*a/3 
        if y <= b/2
            zoneID(e) = 3;
        else 
            zoneID(e) = 6;
        end

    % elseif x <= 4*Lx/6 
    %     if y <= Ly/2
    %         zoneID(e) = 4;
    %     else 
    %         zoneID(e) = 10;
    %     end
    % elseif x <= 5*Lx/6 
    %     if y <= Ly/2
    %         zoneID(e) = 5;
    %     else 
    %         zoneID(e) = 11;
    %     end
    % elseif x <= 6*Lx/6 
    %     if y <= Ly/2
    %         zoneID(e) = 6;
    %     else 
    %         zoneID(e) = 12;
    %     end
     end
end







zone_center = zeros(Nz,2);
for z = 1:Nz
    elems = find(zoneID==z);
    zone_center(z,:) = mean(elemCenter(elems,:),1);
end

Kz = cell(Nz,1);
for z = 1:Nz
    Kz{z} = sparse(Ndof, Ndof);
end
for e = 1:Ne
    z = zoneID(e);
    dofs = elemDof{e};
    Ke = Ke_store{e};
    Kz{z}(dofs,dofs) = Kz{z}(dofs,dofs) + Ke;
end

% check zoning assembly
Ksum = sparse(Ndof,Ndof);
for z=1:Nz, Ksum = Ksum + Kz{z}; end
err = norm(Ksum - K, 'fro') / norm(K, 'fro');
fprintf('Zoning error = %.2e\n', err);




%% assigning the BC for SSSS
Kg=K(Fdof,Fdof);
Mg=M(Fdof,Fdof);
Fg = F(Fdof,1);

% appylying bc to reduced of zonal k
Kzf = cell(Nz,1);
for z = 1:Nz
    Kzf{z} = Kz{z}(Fdof, Fdof);
end

%% Eigen analysis 

[EE,W2]=eigs(Kg,Mg,6,'sm');
omegan=diag(W2).^.5;
frequencies = omegan / (2*pi);  % Convert to Hz

fprintf('\nNatural Frequencies (Hz):\n');
for i = 1:4
    fprintf('Mode %d: %.2f Hz\n', i, frequencies(i));
end


% Rayleigh Damping 
zeta_target = 0.02;
w1 = omegan(1); w2 = omegan(2);
A = [1/w1, w1; 1/w2, w2]; b = 2*zeta_target*[1;1];
alpha_beta = A\b; alphaR = alpha_beta(1); betaR = alpha_beta(2);
Cg = alphaR * Mg + betaR * Kg;

%% -------- model reduction -----------
modes = 4;
Phi_modal = EE(:,1:modes);

% modal reduced matrices
K_modal_red = Phi_modal' * Kg * Phi_modal;
M_modal_red = Phi_modal' * Mg * Phi_modal;
C_modal_red = Phi_modal' * Cg * Phi_modal;
F_modal_red = Phi_modal' * Fg;

% zone stiffness in modal space
Kz_modal_red = cell(Nz,1);
for z = 1:Nz
    Kz_modal_red{z} = Phi_modal' * Kzf{z} * Phi_modal;
end

%% ==================   putting sensor =====================
sensorNode = zeros(Nz,1);
Ns = Nz; % for now 
for z = 1:Nz
    dist = vecnorm(coords - zone_center(z,:), 2, 2);
    [~, sensorNode(z)] = min(dist);
end

% show mesh + zones + sensors
figure('Name','Zones & Sensors'); hold on; axis equal;
cmap = lines(Nz);
for e = 1:Ne
    nodes = conn_nodes(e,:);
    X = coords(nodes,1);
    Y = coords(nodes,2);
    z = zoneID(e);
    patch(X, Y, cmap(z,:), 'EdgeColor','k');
end
for z = 1:Nz
    elems = find(zoneID==z);
    nodes = unique(conn_nodes(elems,:));
    xc = mean(coords(nodes,1)); yc = mean(coords(nodes,2));
    text(xc, yc, num2str(z),'FontWeight','bold');
end
plot(coords(sensorNode,1), coords(sensorNode,2), 'rp','MarkerSize',12,'MarkerFaceColor','r');
title('Zones and Sensors'); xlabel('x'); ylabel('y');
drawnow;


%% ================== NEWMARK-BETA TIME MARCHING =====================
fprintf('\nStarting Newmark-Beta Time Integration...\n');

% 1. Time Simulation Parameters
dt = 0.001;               % Time step (seconds)
T_total = 2.0;             % Total simulation duration (seconds)
time = 0:dt:T_total;       % Time vector
nt = length(time)-1;

%% Apply load function as you want 
load_time_func = sin(2*pi*linspace(0.2,4, nt+1)); 
% Apply load only for the first 0.1s (Impulse-like) if desired, or constant:
% load_time_func(time > 0.1) = 0; 
%% Acting progressive damage factor 
% assign damage factor for different zones 
fd = zeros(Ns,nt+1);
fd(4,:) = linspace(0,0.4,nt+1);
fd(5,:) = linspace(0,0.6,nt+1);
% 3. Newmark-Beta Constants (Average Acceleration Method: Unconditionally Stable)
gamma = 0.5;
beta = 0.25;

% Integration Constants
a0 = 1 / (beta * dt^2);
a1 = gamma / (beta * dt);
a2 = 1 / (beta * dt);
a3 = 1 / (2 * beta) - 1;
a4 = (gamma / beta) - 1;
a5 = (dt / 2) * (gamma / beta - 2);
a6 = dt * (1 - gamma);
a7 = gamma * dt;

% 4. Effective Stiffness Matrix
% K_hat = K + a0*M + a1*C
% Put the initial stiffness matrix
% Kred has the damage term along with it 
Kred_ = sparse(length(Fdof),length(Fdof));
for z=1:Nz
    Kred_ = Kred_ + (1 - fd(z,1))*Kzf{z};
end
K_hat = Kred_ + a0 * Mg + a1 * Cg;

% 5. Initial Conditions (Displacement, Velocity, Acceleration)
% Assume system starts at rest
u_t = zeros(length(Fdof), 1);      % Displacement at current step
ud_t = zeros(length(Fdof), 1);     % Velocity at current step
udd_t = zeros(length(Fdof), 1);    % Acceleration at current step

udd_t = Mg \(Fg*load_time_func(1) - Kred_*u_t -Cg*ud_t );
% Storage for results (Store only sensor node displacement to save memory)
% Find the DOF index for vertical displacement (w) of the sensor node
% The sensorNode index is Global. We need its index in the Reduced system (Fdof).
w_dofs = dof_per_node*((1:Nnode)-1) + 1;
w_dofs_free = intersect(w_dofs,Fdof);
sensor_global_dof = (sensorNode-1)*dof_per_node + 1; % DOF 1 is 'w'
[found, sensor_reduced_idx] = ismember(sensor_global_dof, w_dofs_free);
    
    %sensor_reduced_idx = find(Fdof == sensor_global_dof); % this the what we get after reducing the constrained dofs 

sensor_disp_history_true = zeros(Ns, nt+1);
sensor_vel_history_true = zeros(Ns, nt+1);
sensor_accn_hist_true = zeros(Ns,nt+1);
sensor_accn_hist_true(:,1) = udd_t(sensor_reduced_idx);
% 6. Time Stepping Loop
for i = 2:nt+1
    % Current Load Vector at time t+dt
    % F_dynamic = F_static_shape * time_function
    F_t_plus_dt = Fg * load_time_func(i);

    % Calculating the changing stiffness with time 
    Kred = sparse(length(Fdof),length(Fdof));
    for z=1:Nz
        Kred = Kred + (1 - fd(z,i))*Kzf{z};
    end

    % Calculate Effective Load at t+dt
    % R_hat = F(t+dt) + M*(a0*u + a2*ud + a3*udd) + C*(a1*u + a4*ud + a5*udd)
    
    M_term = Mg * (a0 * u_t + a2 * ud_t + a3 * udd_t);
    C_term = Cg * (a1 * u_t + a4 * ud_t + a5 * udd_t);
    
    R_hat = F_t_plus_dt + M_term + C_term;
    
    % Solve for Displacement at t+dt
    % K_hat * u(t+dt) = R_hat
    K_hat = Kred + a0 * Mg + a1 * Cg;
    u_next = K_hat \ R_hat;
    
    % Update Velocity and Acceleration at t+dt
    udd_next = a0 * (u_next - u_t) - a2 * ud_t - a3 * udd_t;
    ud_next = ud_t + a6 * udd_t + a7 * udd_next;
    
    % Store Sensor Data
    if ~isempty(sensor_reduced_idx)
        sensor_disp_history_true(:,i) = u_next(sensor_reduced_idx);
        sensor_accn_hist_true(:,i) = udd_next(sensor_reduced_idx);
        sensor_vel_history_true(:,i) = ud_next(sensor_reduced_idx);
    end
    
    % Update states for next iteration
    u_t = u_next;
    ud_t = ud_next;
    udd_t = udd_next;
end

fprintf('Time integration complete.\n');
 
%% ---- Build mapping Phi_w and sensor selection S ----
% rows in Phi_modal correspond to free DOFs only (size nFree x modes)
% we need rows corresponding to free w DOFs
w_dof = 4*((1:Nnode)-1) +1;
w_is = ismember(Fdof, w_dof);        % logical index into freeDOF for all w DOFs
% Phi_w = Phi_modal(w_is, :);              % Nwfree x modes

% % sensor DOFs: global w DOFs of sensor nodes
% sensor_w_global = 3*sensorNode - 2;      % size Nz x 1 (global DOF index)
% % convert to index among free w DOFs
% [found, sensor_idx] = ismember(sensor_w_global, w_dofs_free);
% if ~all(found)
%     error('One or more sensor nodes are constrained (on boundary). Move sensor.');
% end

Ns = length(sensor_reduced_idx);
S = zeros(Ns, length(Fdof));
for i = 1:Ns
    S(i, sensor_reduced_idx(i)) = 1;
end

% sensor-mode matrix (maps modal coords to sensor w)
Phi_sen = S * Phi_modal;   % size Ns x modes


%% ------- Now the Time for DEKF ----------------------
params = Ns;
x_state = zeros(2*modes,1); % initial state
fd_param = fd(:,1); % initial damage factor 

% ADDING NOISE to measured value
noise_a = 1e-8;
noise_u = 1e-8;
noise_v = 1e-8;

% Initialize the state covariance matrix and process noise covariance
Px = 1e-6*eye(2*modes); % State covariance matrix
Qx = 5e-6 * eye(2*modes); % Process noise covariance
Pfd = 1e-11*eye(params);
Qfd = 1e-11*eye(params);

R_u = 1e-9 * eye(params);
R_v = 1e-9 * eye(params);
Rx = blkdiag(R_u, R_v);           % measurement covariance for u and v (2*Ns x 2*Ns)
Rfd = 1e-2 * eye(params); 

% store the solution for future use 
u_dkf_hist = zeros(Ns,nt+1);
v_dkf_hist = zeros(Ns,nt+1);  % contains only the sensor's states 
a_dkf_hist = zeros(Ns, nt+1);
a_dkf_hist_state = zeros(Ns, nt+1);

x_state_hist = zeros(2*modes,nt+1);
fd_est_hist = zeros(params,nt+1);
% Initialize state history
x_state_hist(:, 1) = x_state;
fd_est_hist(:, 1) = fd_param;


% starting the loop 
for i=1:nt
     F_t = F_modal_red * load_time_func(i);
    % getting the measurement accn 
    a_meas = sensor_accn_hist_true(:,i) + noise_a*randn;
    u_meas = sensor_disp_history_true(:,i) + noise_u*randn;
    v_meas = sensor_vel_history_true(:,i) + noise_v * randn;

     K_dkf_modal = zeros(modes,modes);
    for z = 1:Nz
        K_dkf_modal = K_dkf_modal + (1 - fd_param(z)) * Kz_modal_red{z};
    end


     a_dkf_prev =  M_modal_red \ (F_t - K_dkf_modal*x_state(1:modes) - C_modal_red*x_state(modes+1:end));

     % random walk prediction
     fd_pred = fd_param + 1e-8*diag(params);
     Pfd_pred = Pfd + Qfd;
    
     %now using predicted Fd
      K_dkf_modal = zeros(modes,modes);
    for z = 1:Nz
        K_dkf_modal = K_dkf_modal + (1 - fd_pred(z)) * Kz_modal_red{z};
    end


     H_fd = zeros(Ns, params);
    for z = 1:params
        % derivative: S * M^{-1} * Kz * u
        % using modal: Phi_sen * (M_modal_red \ (Kz_modal_red{z} * q_pred))
        tmp = (M_modal_red + dt*0.5*C_modal_red + dt^2*0.25*K_dkf_modal)\ (Kz_modal_red{z} * (x_state(1:modes) + dt*x_state(modes+1:end) + dt^2*0.25*a_dkf_prev));    % modes x 1
        H_fd(:, z) = Phi_sen * tmp;                        % Ns x 1
    end
    % using the Eulerian scheme
    

    a_dkf_pred = (M_modal_red + dt*0.5*C_modal_red + dt^2*0.25*K_dkf_modal)\(F_modal_red*load_time_func(i+1) - K_dkf_modal*(x_state(1:modes) + dt*x_state(modes+1:end) + dt^2*0.25*a_dkf_prev) - C_modal_red* (x_state(modes+1:end) + dt*0.5*a_dkf_prev)  ); 
    a_dkf_hist(:,i+1) = Phi_sen*a_dkf_pred;
    % Updation part 
    z_accn_diff = a_meas - Phi_sen*a_dkf_pred;

      % update parameters using EKF formula
    S_fd = H_fd * Pfd_pred * H_fd' + Rfd;
    S_fd = S_fd + 1e-12 * eye(size(S_fd));  % regularize
    K_fd = (Pfd_pred * H_fd') / S_fd;      % Nz x Ns

    fd_upd = fd_pred + K_fd * z_accn_diff;
    % enforce physical bounds and remove NaNs
    fd_upd = real(fd_upd);                  % remove tiny imaginary parts
    fd_upd(isnan(fd_upd)) = 0;              % safety
    fd_upd = max(0, min(0.99, fd_upd));     % clamp to [0,0.99]

    Pfd = (eye(params) - K_fd * H_fd) * Pfd_pred;

    fd_param = fd_upd;
    fd_est_hist(:, i+1) = fd_param;

    %% Now for coupling with the states 
    % state will be updated for the next step (i+1) 
    Ft_plus = F_modal_red*load_time_func(i+1); 
     K_dkf_modal = zeros(modes,modes);
    for z = 1:Nz
        K_dkf_modal = K_dkf_modal + (1 - fd_param(z)) * Kz_modal_red{z};
    end

     % Continuous-time modal A_cont (correct sign)
    A_cont = [ zeros(modes), eye(modes);
              - (M_modal_red \ K_dkf_modal),  - (M_modal_red \ C_modal_red) ];  % 2m x 2m

    % Continuous-time B_cont maps modal force (F_modal_red) into states
    B_cont = [ zeros(modes); (M_modal_red \ eye(modes)) ];   % 2m x m
    % Discretize (exact) using matrix exponential on augmented matrix
    nA = 2*modes;
    Mbig = [A_cont, B_cont; zeros(modes, nA), zeros(modes)];
    E = expm(Mbig*dt);
    Ad = E(1:nA, 1:nA);
    Bd = E(1:nA, nA+1:end);

    % predict
    x_pred = Ad * x_state + Bd * Ft_plus;
    Px_pred = Ad * Px * Ad' + Qx;

    % updation of the states 
    Hx = [Phi_sen, zeros(Ns, modes); zeros(Ns, modes), Phi_sen];  % 2Ns x 2m

    Z_uv = [u_meas; v_meas];    % 2Ns x 1
    z_uv_diff = Z_uv - Hx * x_pred;

    Sx = Hx * Px_pred * Hx' + Rx;
    % regularize Sx to avoid singularity
    Sx = Sx + 1e-12 * eye(size(Sx));
    Kx = (Px_pred * Hx') / Sx;       % 2m x 2Ns

    x_upd = x_pred + Kx * z_uv_diff;
    Px = (eye(2*modes) - Kx * Hx) * Px_pred;

    x_state = x_upd;
    x_state_hist(:, i+1) = x_state;
    u_dkf_hist(:,i+1) = Phi_sen*x_state(1:modes);

    a_dkf_hist_state(:,i+1) = Phi_sen* (M_modal_red \ (Ft_plus - K_dkf_modal*x_state(1:modes) - C_modal_red*x_state(modes+1:end)) ) ;   

end






%% ================== PLOTTING RESPONSE =====================
figure('Name', 'Time History Response');
for i=1:Ns
    subplot(2,3,i);
    plot(time, sensor_disp_history_true(i,:), 'k-', 'LineWidth', 1.5); hold on;
    plot(time, u_dkf_hist(i,:), 'b--', 'LineWidth', 1.5);
    grid on;
    title(['Vertical Displacement at Sensor Node ' num2str(sensorNode(i))]); 
    xlabel('Time (s)');
    ylabel('Displacement (m)');
    xlim([0, T_total]);
end


% for accn
figure('Name', 'Time History Response');
for i=1:Ns
    subplot(2,3,i);
    plot(time, sensor_accn_hist_true(i,:), 'k-', 'LineWidth', 1.5); hold on;
    plot(time, a_dkf_hist(i,:) ,'b--', 'LineWidth', 1.5);
    plot(time, a_dkf_hist_state(i,:), 'm--', 'LineWidth', 1.5);
    grid on;
    title(['Vertical Accn at Sensor Node ' num2str(sensorNode(i))]);
    xlabel('Time (s)');
    ylabel('Accn (m/s2)');
    xlim([0, T_total]);
end


% for damage detection plot
figure('Name', ' Prgressive Damage Response');
for i=1:Ns
    subplot(2,3,i);
    plot(time, fd(i,:), 'k-', 'LineWidth', 1.5); hold on;
    plot(time, fd_est_hist(i,:), 'b--', 'LineWidth', 1.5);
    grid on;
    title(['Damage detection  ' num2str(sensorNode(i))]);
    xlabel('Time (s)');
    ylabel('damage Factor');
    xlim([0, T_total]);
end

