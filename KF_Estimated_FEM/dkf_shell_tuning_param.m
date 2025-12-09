%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Souvik Sarkar 
 %                22ME10083
 %    DKF _shell but needd to tune P,Q, R for optimaum result
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc;

%% ------------------- MATERIAL & GEOMETRY ------------------------------
E = 210e9;                      % Young's modulus (Pa)
t = 0.25;                       % Thickness (m)
rho = 8050;                     % density (kg/m3)
q = -2e6;                       % Transverse load (N/m2) (example)
nu = 0.3;
mu = 5/6;
global_z = [0 0 1];

D = (E / (1 - nu^2)) * [1, nu, 0, 0, 0, 0;
                        nu, 1, 0, 0, 0, 0;
                        0, 0, 0, 0, 0, 0;
                        0, 0, 0, (1 - nu)/2, 0, 0;
                        0, 0, 0, 0,  mu*((1 - nu)/2), 0;
                        0, 0, 0, 0, 0, mu*((1 - nu)/2)];

%% ------------------- NODES & CONNECTIVITY -----------------------------
coords = [ ...
 0,0,0, 0,0,0.2;
 0.25,0,0, 0.25,0,0.2;
 0.5,0,0, 0.5,0,0.2;
 0.75,0,0, 0.75,0,0.2;
 1,0,0, 1,0,0.2;
 1.25,0,0, 1.25,0,0.2;
 1.5,0,0, 1.5,0,0.2;
 1.75,0,0, 1.75,0,0.2;
 2,0,0, 2,0,0.2;

 0,0.25,0, 0,0.25,0.2;
 0.5,0.25,0, 0.5,0.25,0.2;
 1,0.25,0, 1,0.25,0.2;
 1.5,0.25,0, 1.5,0.25,0.2;
 2,0.25,0, 2,0.25,0.2;

 0,0.5,0, 0,0.5,0.2;
 0.25,0.5,0, 0.25,0.5,0.2;
 0.5,0.5,0, 0.5,0.5,0.2;
 0.75,0.5,0, 0.75,0.5,0.2;
 1,0.5,0, 1,0.5,0.2;
 1.25,0.5,0, 1.25,0.5,0.2;
 1.5,0.5,0, 1.5,0.5,0.2;
 1.75,0.5,0, 1.75,0.5,0.2;
 2,0.5,0, 2,0.5,0.2;

 0,0.75,0, 0,0.75,0.2;
 0.5,0.75,0, 0.5,0.75,0.2;
 1,0.75,0, 1,0.75,0.2;
 1.5,0.75,0, 1.5,0.75,0.2;
 2,0.75,0, 2,0.75,0.2;

 0,1,0, 0,1,0.2;
 0.25,1,0, 0.25,1,0.2;
 0.5,1,0, 0.5,1,0.2;
 0.75,1,0, 0.75,1,0.2;
 1,1,0, 1,1,0.2;
 1.25,1,0, 1.25,1,0.2;
 1.5,1,0, 1.5,1,0.2;
 1.75,1,0, 1.75,1,0.2;
 2,1,0, 2,1,0.2;

 0,1.25,0, 0,1.25,0.2;
 0.5,1.25,0, 0.5,1.25,0.2;
 1,1.25,0, 1,1.25,0.2;
 1.5,1.25,0, 1.5,1.25,0.2;
 2,1.25,0, 2,1.25,0.2;

 0,1.5,0, 0,1.5,0.2;
 0.25,1.5,0, 0.25,1.5,0.2;
 0.5,1.5,0, 0.5,1.5,0.2;
 0.75,1.5,0, 0.75,1.5,0.2;
 1,1.5,0, 1,1.5,0.2;
 1.25,1.5,0, 1.25,1.5,0.2;
 1.5,1.5,0, 1.5,1.5,0.2;
 1.75,1.5,0, 1.75,1.5,0.2;
 2,1.5,0, 2,1.5,0.2;

 0,1.75,0, 0,1.75,0.2;
 0.5,1.75,0, 0.5,1.75,0.2;
 1,1.75,0, 1,1.75,0.2;
 1.5,1.75,0, 1.5,1.75,0.2;
 2,1.75,0, 2,1.75,0.2;

 0,2,0, 0,2,0.2;
 0.25,2,0, 0.25,2,0.2;
 0.5,2,0, 0.5,2,0.2;
 0.75,2,0, 0.75,2,0.2;
 1,2,0, 1,2,0.2;
 1.25,2,0, 1.25,2,0.2;
 1.5,2,0, 1.5,2,0.2;
 1.75,2,0, 1.75,2,0.2;
 2,2,0, 2,2,0.2
];

s = size(coords);

c = [ ...
 1,  3, 17, 15,  2, 11, 16, 10;
 3,  5, 19, 17,  4, 12, 18, 11;
 5,  7, 21, 19,  6, 13, 20, 12;
 7,  9, 23, 21,  8, 14, 22, 13;

15, 17, 31, 29, 16, 25,30, 24;
17, 19, 33, 31, 18, 26, 32, 25;
19, 21, 35, 33, 20, 27, 34, 26;
21, 23, 37,35, 22, 28, 36, 27;

29, 31, 45, 43, 30, 39, 44, 38;
31, 33, 47, 45, 32, 40, 46, 39;
33, 35, 49, 47, 34, 41, 48, 40;
35, 37, 51, 49, 36, 42, 50, 41;

43, 45, 59, 57, 44, 53, 58, 52;
45, 47, 61, 59, 46, 54, 60, 53;
47, 49, 63, 61, 48, 55, 62, 54;
49, 51, 65, 63, 50, 56, 64, 55
];

Nelem = 16;
Nnode = s(1);
Ndof = 5*Nnode;

%% midpoints
mid_coords=zeros(Nnode,3);
for i = 1:Nnode
    mid_coords(i,:) = 0.5*[coords(i,1)+coords(i,4), coords(i,2)+coords(i,5), coords(i,3)+coords(i,6)];
end

%% shape functions (symbolic once)
syms xi eta real
N1 = -0.25*(1 - xi)*(1 - eta)*(1 + xi + eta);
N2 = -0.25*(1 + xi)*(1 - eta)*(1 - xi + eta);
N3 = -0.25*(1 + xi)*(1 + eta)*(1 - xi - eta);
N4 = -0.25*(1 - xi)*(1 + eta)*(1 + xi - eta);
N5 =  0.5*(1 - xi^2)*(1 - eta);
N6 =  0.5*(1 + xi)*(1 - eta^2);
N7 =  0.5*(1 - xi^2)*(1 + eta);
N8 =  0.5*(1 - xi)*(1-eta^2);
N = [N1,N2,N3,N4,N5,N6,N7,N8];

N_xi = simplify(diff(N, xi));
N_eta = simplify(diff(N, eta));
dN_dxi = [N_xi; N_eta];

%% Gauss points & weights (8-point)
gp = [-sqrt(0.6), -sqrt(0.6);
       sqrt(0.6), -sqrt(0.6);
       sqrt(0.6),  sqrt(0.6);
      -sqrt(0.6),  sqrt(0.6);
       0, -sqrt(0.6);
       sqrt(0.6), 0;
       0, sqrt(0.6);
      -sqrt(0.6), 0];
w = [25/81,25/81,25/81,25/81,40/81,40/81,40/81,64/81]';

%% Global assemble
Kg=zeros(Ndof,Ndof);
Kg_undam = zeros(Ndof,Ndof);
Mg=zeros(Ndof,Ndof);
Fg = zeros(Ndof,1);

nodeXiEta = [-1,-1; 0,-1; 1,-1; 1,0; 1,1; 0,1; -1,1; -1,0];

theta_ori = ones(Nelem,1);  % initial (no damage)

fprintf('Assembling element matrices...\n');
for e = 1:Nelem
    nodes = c(e,:);
    XYZ = mid_coords(nodes,:);
    Nof = zeros(3,40); Ntf=zeros(3,40);
    Bzetaf=zeros(6,40); Bof=zeros(6,40);
    Ke=zeros(40,40); Me=zeros(40,40); Fe=zeros(40,1);
    V_edge=zeros(8,3);
    for cnt=1:8
        k = nodes(cnt);
        V_edge(cnt,:) = -[(coords(k,1)-coords(k,4)), (coords(k,2)-coords(k,5)), (coords(k,3)-coords(k,6))];
    end
    for i = 1:8
        xi_val = gp(i,1); eta_val = gp(i,2); wts = w(i);
        N_val = double(subs(N, [xi, eta], [xi_val, eta_val]));
        dN_dxi_val = double(subs(dN_dxi, [xi, eta], [xi_val, eta_val]));
        dx_dxi = dN_dxi_val * XYZ;
        J = zeros(3,3); J(1:2,:) = dx_dxi; J(3,:) = 0.5 * N_val * V_edge;
        J_ac = J'; detJ_ac = det(J_ac);
        if abs(detJ_ac) < 1e-12, detJ_ac = sign(detJ_ac)*1e-12; end
        J_ = inv(J');  % okay because detJ handled
        Jac = [J_(1,1), J_(2,1), J_(3,1), 0,0,0,0,0,0;
               0,0,0, J_(1,2), J_(2,2), J_(3,2),0,0,0;
               0,0,0,0,0,0, J_(1,3), J_(2,3), J_(3,3);
               J_(1,2), J_(2,2), J_(3,2), J_(1,1), J_(2,1), J_(3,1),0,0,0;
               0,0,0, J_(1,3), J_(2,3), J_(3,3), J_(1,2), J_(2,2), J_(3,2);
               J_(1,3), J_(2,3), J_(3,3),0,0,0, J_(1,1), J_(2,1), J_(3,1)];
        for cnt = 1:8
            k = nodes(cnt);
            xi_ = nodeXiEta(cnt,1); eta_ = nodeXiEta(cnt,2);
            dN_dxi_ = double(subs(dN_dxi, [xi, eta], [xi_, eta_]));
            e1 = dN_dxi_(1,:) * XYZ; e2 = dN_dxi_(2,:) * XYZ;
            e1_unit = e1 / norm(e1); e2_unit = e2 / norm(e2);
            e3 = cross(e1_unit, e2_unit); e3_unit = e3 / norm(e3);
            v3 = e3_unit; v1 = e1_unit;
            if dot(v3, global_z) < 0
                v3 = -v3;
            end
            v2 = cross(v3, v1); v2 = v2 / norm(v2);
            lx = v1(1); ly = v2(1); lz = v3(1);
            mx = v1(2); my = v2(2); mz = v3(2);
            nx = v1(3); ny = v2(3); nz = v3(3);
            T_eps = [lx^2, mx^2, nx^2, 2*lx*mx, 2*mx*nx, 2*nx*lx;
                     ly^2, my^2, ny^2, 2*ly*my, 2*my*ny, 2*ny*ly;
                     lz^2, mz^2, nz^2, 2*lz*mz, 2*mz*nz, 2*nz*lz;
                     lx*ly, mx*my, nx*ny, lx*my+mx*ly, mx*ny+nx*my, nx*ly+lx*ny;
                     ly*lz, my*mz, ny*nz, ly*mz+my*lz, my*nz+ny*mz, ny*lz+ly*nz;
                     lz*lx, mz*mx, nz*nx, lz*mx+mz*lx, mz*nx+nz*mx, nz*lx+lz*nx];
            Nof(:,(5*(cnt-1)+1):5*cnt)= N_val(cnt)*[1,0,0,0,0; 0,1,0,0,0; 0,0,1,0,0];
            Ntf(:,(5*(cnt-1)+1):5*cnt)= N_val(cnt)*[0,0,0,-0.5*t*v2(1),0.5*t*v1(1);
                                                  0,0,0,-0.5*t*v2(2),0.5*t*v1(2);
                                                  0,0,0,-0.5*t*v2(3),0.5*t*v1(3)];
            E_mat = T_eps * Jac;
            % Bof using the quadrature evaluated dN_dxi_val
            Bof(:,(5*(cnt-1)+1):5*cnt) = E_mat * [dN_dxi_val(1,cnt), 0, 0, 0, 0;
                                                   dN_dxi_val(2,cnt), 0, 0, 0, 0;
                                                   0, 0, 0, -N_val(cnt)*0.5*t*v2(1), N_val(cnt)*0.5*t*v1(1);
                                                   0, dN_dxi_val(1,cnt), 0, 0, 0;
                                                   0, dN_dxi_val(2,cnt), 0, 0, 0;
                                                   0, 0, 0, -N_val(cnt)*0.5*t*v2(2), N_val(cnt)*0.5*t*v1(2);
                                                   0, 0, dN_dxi_val(1,cnt), 0, 0;
                                                   0, 0, dN_dxi_val(2,cnt), 0, 0;
                                                   0, 0, 0, -N_val(cnt)*0.5*t*v2(3), N_val(cnt)*0.5*t*v1(3)];
            Bzetaf(:,(5*(cnt-1)+1):5*cnt) = E_mat * [0,0,0,-dN_dxi_val(1,cnt)*t*0.5*v2(1),dN_dxi_val(1,cnt)*t*0.5*v1(1);
                                                     0,0,0,-dN_dxi_val(2,cnt)*t*0.5*v2(1),dN_dxi_val(2,cnt)*t*0.5*v1(1);
                                                     0,0,0,0,0;
                                                     0,0,0,-dN_dxi_val(1,cnt)*t*0.5*v2(2),dN_dxi_val(1,cnt)*t*0.5*v1(2);
                                                     0,0,0,-dN_dxi_val(2,cnt)*t*0.5*v2(2),dN_dxi_val(2,cnt)*t*0.5*v1(2);
                                                     0,0,0,0,0;
                                                     0,0,0,-dN_dxi_val(1,cnt)*t*0.5*v2(3),dN_dxi_val(1,cnt)*t*0.5*v1(3);
                                                     0,0,0,-dN_dxi_val(2,cnt)*t*0.5*v2(3),dN_dxi_val(2,cnt)*t*0.5*v1(3);
                                                     0,0,0,0,0];
            Fe(5*(cnt-1)+3) = Fe(5*(cnt-1)+3) + q * N_val(cnt) * abs(detJ_ac) * wts;
        end
        Ke = Ke + ((2*Bof'*D*Bof) + (2/3)*Bzetaf'*D*Bzetaf) * abs(detJ_ac) * wts;
        Me = Me + rho * (2*(Nof)'*Nof + (2/3)*(Ntf')*Ntf) * abs(detJ_ac) * wts;
    end
    dofs = reshape([5*(nodes-1)+1; 5*(nodes-1)+2; 5*(nodes-1)+3; 5*(nodes-1)+4; 5*nodes], 1, []);
    Kg(dofs,dofs) = Kg(dofs,dofs) + theta_ori(e)*Ke;
    Kg_undam(dofs,dofs) = Kg_undam(dofs,dofs) + Ke;
    Mg(dofs,dofs) = Mg(dofs,dofs) + Me;
    Fg(dofs,1) = Fg(dofs,1) + Fe;
end
fprintf('Assembly done.\n');

%% ------------------ Boundary conditions & reduce -----------------------
fixedDOF = [1,2,3,4,5,46,47,48,49,50,71,72,73,74,75,116,117,118,119,120,141,142,143,144,145,...
            186,187,188,189,190,211,212,213,214,215,256,257,258,259,260,281,282,283,284,285];
freeDOF = setdiff(1:Ndof, fixedDOF);

Kf = Kg(freeDOF,freeDOF);
Kf_undam = Kg_undam(freeDOF,freeDOF);
Mf = Mg(freeDOF,freeDOF);
Ff = Fg(freeDOF);

%% ----------------------- Static check & plotting -----------------------
U = zeros(Ndof,1);
U(freeDOF) = Kf \ Ff;
U_undam = zeros(Ndof,1);
U_undam(freeDOF) = Kf_undam \ Ff;

%% -------------------- Modal analysis & Rayleigh -------------------------
nEig = 10;
try
    [PHI_all, W2_all] = eigs(Kf_undam, Mf, nEig, 'smallestabs');
    omega_all = sqrt(abs(diag(W2_all)));
catch
    [PHIs, W2s] = eig(Kf_undam, Mf);
    omega_all = sqrt(abs(diag(W2s)));
    [omega_all, idxs] = sort(omega_all);
    PHI_all = PHIs(:, idxs(1:nEig));
    omega_all = omega_all(1:nEig);
end
freq = omega_all / (2*pi);
zeta_target = 0.02;
w1 = omega_all(1); w2 = omega_all(min(3,length(omega_all)));
A = [1/w1, w1; 1/w2, w2]; b = 2*zeta_target * [1;1];
alpha_beta = A\b; alphaR = alpha_beta(1); betaR = alpha_beta(2);
Cf = alphaR * Mf + betaR * Kf_undam;

%% -------------------- Modal reduction for filters ----------------------
n_modes = min(4, length(omega_all));
PHI_red = PHI_all(:,1:n_modes);
for i = 1:n_modes
    PHI_red(:,i) = PHI_red(:,i) / sqrt(PHI_red(:,i)' * Mf * PHI_red(:,i));
end
K_modal = PHI_red' * Kf_undam * PHI_red;
M_modal = PHI_red' * Mf * PHI_red;
C_modal = PHI_red' * Cf * PHI_red;
F_modal = PHI_red' * Ff;

fprintf('Modal reduction: %d modes\n', n_modes);



%% ------------------ Time integration (Newmark) to produce true signals ----
beta = 1/4; gamma = 1/2;
dt = 0.001; Ttotal = 1.0; tvec = 0:dt:Ttotal; nSteps = length(tvec)-1;
a0 = 1/(beta*dt^2); a1 = gamma/(beta*dt); a2 = 1/(beta*dt); a3 = 1/(2*beta)-1; a4 = gamma/beta - 1; a5 = dt*(gamma/(2*beta)-1);
% -------- 
% Generating damage 
theta_true = linspace(0, 0.5, nSteps+1);
%for i



% We'll simulate a time-varying damage by scaling Kf_undam using theta profile
% define element-level theta evolution (time x elements) and assemble K_t each step
% Ntime = nSteps+1;
% theta_true = ones(Nelem, Ntime);
% % Example: elements 2,6,7,8,10 progressively drop from 1.0 -> 0.6 between t0 and t1
% damaged_elems = [2,6,7,8,10];
% t0 = 0.2; t1 = 0.6; theta_min = 0.6;
% for kt = 1:Ntime
%     t = tvec(kt);
%     for e = damaged_elems
%         if t < t0
%             theta_true(e,kt) = 1.0;
%         elseif t > t1
%             theta_true(e,kt) = theta_min;
%         else
%             theta_true(e,kt) = 1.0 - (1.0 - theta_min)*((t - t0)/(t1 - t0));
%         end
%     end
% end

% To apply theta to global stiffness at each time step we need mapping element->global DOFs.
% Re-assemble global stiffness per time-step is expensive, but since Nelem small (16) it's fine.
% Pre-store undamaged element stiffness contributions (Ke per element and dof indices) from assembly above.
% For that, we need to re-run element loop but store Ke and corresponding dof arrays.
fprintf('Extracting elemental Ke/Me/dofs (for runtime theta scaling)...\n');
ElementKe = cell(Nelem,1);
ElementMe = cell(Nelem,1);
ElementDofs = cell(Nelem,1);
% We'll recompute element matrices (similar to earlier code) but now store Ke/Me/dofs
for e = 1:Nelem
    nodes = c(e,:);
    XYZ = mid_coords(nodes,:);
    Nof = zeros(3,40); Ntf=zeros(3,40);
    Bzetaf=zeros(6,40); Bof=zeros(6,40);
    Ke=zeros(40,40); Me=zeros(40,40);
    V_edge=zeros(8,3);
    for cnt=1:8
        k = nodes(cnt);
        V_edge(cnt,:) = -[(coords(k,1)-coords(k,4)), (coords(k,2)-coords(k,5)), (coords(k,3)-coords(k,6))];
    end
    for i = 1:8
        xi_val = gp(i,1); eta_val = gp(i,2); wts = w(i);
        N_val = double(subs(N, [xi, eta], [xi_val, eta_val]));
        dN_dxi_val = double(subs(dN_dxi, [xi, eta], [xi_val, eta_val]));
        dx_dxi = dN_dxi_val * XYZ;
        J = zeros(3,3); J(1:2,:) = dx_dxi; J(3,:) = 0.5 * N_val * V_edge;
        J_ac = J'; detJ_ac = det(J_ac);
        if abs(detJ_ac) < 1e-12, detJ_ac = sign(detJ_ac)*1e-12; end
        J_ = inv(J'); Jac = [J_(1,1), J_(2,1), J_(3,1), 0,0,0,0,0,0;
               0,0,0, J_(1,2), J_(2,2), J_(3,2),0,0,0;
               0,0,0,0,0,0, J_(1,3), J_(2,3), J_(3,3);
               J_(1,2), J_(2,2), J_(3,2), J_(1,1), J_(2,1), J_(3,1),0,0,0;
               0,0,0, J_(1,3), J_(2,3), J_(3,3), J_(1,2), J_(2,2), J_(3,2);
               J_(1,3), J_(2,3), J_(3,3),0,0,0, J_(1,1), J_(2,1), J_(3,1)];
        for cnt = 1:8
            k = nodes(cnt);
            dN_dxi_val_c = dN_dxi_val;
            e1 = double(subs(N_xi, [xi,eta],[nodeXiEta(cnt,1), nodeXiEta(cnt,2)])) * XYZ; %#ok<NASGU>
            % (we re-use earlier expressions; keep consistent)
            % compute local axes:
            xi_ = nodeXiEta(cnt,1); eta_ = nodeXiEta(cnt,2);
            dN_dxi_ = double(subs(dN_dxi, [xi, eta], [xi_, eta_]));
            e1 = dN_dxi_(1,:) * XYZ; e2 = dN_dxi_(2,:) * XYZ;
            e1_unit = e1 / norm(e1); e2_unit = e2 / norm(e2);
            e3 = cross(e1_unit, e2_unit); e3_unit = e3 / norm(e3);
            v3 = e3_unit; v1 = e1_unit;
            if dot(v3,global_z) < 0, v3 = -v3; end
            v2 = cross(v3,v1); v2 = v2 / norm(v2);
            T_eps = zeros(6,6);
            lx = v1(1); ly = v2(1); lz = v3(1);
            mx = v1(2); my = v2(2); mz = v3(2);
            nx = v1(3); ny = v2(3); nz = v3(3);
            T_eps = [lx^2, mx^2, nx^2, 2*lx*mx, 2*mx*nx, 2*nx*lx;
                     ly^2, my^2, ny^2, 2*ly*my, 2*my*ny, 2*ny*ly;
                     lz^2, mz^2, nz^2, 2*lz*mz, 2*mz*nz, 2*nz*lz;
                     lx*ly, mx*my, nx*ny, lx*my+mx*ly, mx*ny+nx*my, nx*ly+lx*ny;
                     ly*lz, my*mz, ny*nz, ly*mz+my*lz, my*nz+ny*mz, ny*lz+ly*nz;
                     lz*lx, mz*mx, nz*nx, lz*mx+mz*lx, mz*nx+nz*mx, nz*lx+lz*nx];
            Nof(:,(5*(cnt-1)+1):5*cnt)= N_val(cnt)*[1,0,0,0,0; 0,1,0,0,0; 0,0,1,0,0];
            Ntf(:,(5*(cnt-1)+1):5*cnt)= N_val(cnt)*[0,0,0,-0.5*t*v2(1),0.5*t*v1(1);
                                                  0,0,0,-0.5*t*v2(2),0.5*t*v1(2);
                                                  0,0,0,-0.5*t*v2(3),0.5*t*v1(3)];
            E_mat = T_eps * Jac;
            Bof(:,(5*(cnt-1)+1):5*cnt) = E_mat * [dN_dxi_val(1,cnt), 0, 0, 0, 0;
                                                  dN_dxi_val(2,cnt), 0, 0, 0, 0;
                                                  0, 0, 0, -N_val(cnt)*0.5*t*v2(1), N_val(cnt)*0.5*t*v1(1);
                                                  0, dN_dxi_val(1,cnt), 0, 0, 0;
                                                  0, dN_dxi_val(2,cnt), 0, 0, 0;
                                                  0, 0, 0, -N_val(cnt)*0.5*t*v2(2), N_val(cnt)*0.5*t*v1(2);
                                                  0, 0, dN_dxi_val(1,cnt), 0, 0;
                                                  0, 0, dN_dxi_val(2,cnt), 0, 0;
                                                  0, 0, 0, -N_val(cnt)*0.5*t*v2(3), N_val(cnt)*0.5*t*v1(3)];
            Bzetaf(:,(5*(cnt-1)+1):5*cnt) = E_mat * [0,0,0,-dN_dxi_val(1,cnt)*t*0.5*v2(1),dN_dxi_val(1,cnt)*t*0.5*v1(1);
                                                     0,0,0,-dN_dxi_val(2,cnt)*t*0.5*v2(1),dN_dxi_val(2,cnt)*t*0.5*v1(1);
                                                     0,0,0,0,0;
                                                     0,0,0,-dN_dxi_val(1,cnt)*t*0.5*v2(2),dN_dxi_val(1,cnt)*t*0.5*v1(2);
                                                     0,0,0,-dN_dxi_val(2,cnt)*t*0.5*v2(2),dN_dxi_val(2,cnt)*t*0.5*v1(2);
                                                     0,0,0,0,0;
                                                     0,0,0,-dN_dxi_val(1,cnt)*t*0.5*v2(3),dN_dxi_val(1,cnt)*t*0.5*v1(3);
                                                     0,0,0,-dN_dxi_val(2,cnt)*t*0.5*v2(3),dN_dxi_val(2,cnt)*t*0.5*v1(3);
                                                     0,0,0,0,0];
            % accumulate Ke/Me for this gauss pt
        end
        Ke = Ke + ((2*Bof'*D*Bof) + (2/3)*Bzetaf'*D*Bzetaf) * abs(detJ_ac) * wts;
        Me = Me + rho * (2*(Nof)'*Nof + (2/3)*(Ntf')*Ntf) * abs(detJ_ac) * wts;
    end
    dofs = reshape([5*(nodes-1)+1; 5*(nodes-1)+2; 5*(nodes-1)+3; 5*(nodes-1)+4; 5*nodes], 1, []);
    ElementKe{e} = Ke;
    ElementMe{e} = Me;
    ElementDofs{e} = dofs;
end
fprintf('Element Ke extraction done.\n');

% Newmark integration with time-varying theta applied to K
Kf_base = Kf_undam; % base undamaged reduced-global stiffness (free DOF)
nFree = size(Kf_base,1);
u_true = zeros(nFree,1); v_true=zeros(nFree,1);
a_true = Mf \ (Ff * 0 - Cf*v_true - Kf_base*u_true); % initial acceleration

Ufree_hist_true = zeros(nFree, nSteps+1);
Vfree_hist_true = zeros(nFree, nSteps+1);
Afree_hist_true = zeros(nFree, nSteps+1);
Utotal_hist_true = zeros(5*Nnode, nSteps+1);

Ufree_hist_true(:,1) = u_true; Vfree_hist_true(:,1) = v_true; Afree_hist_true(:,1) = a_true;
Utemp = zeros(5*Nnode,1); Utemp(freeDOF) = u_true; Utotal_hist_true(:,1) = Utemp;

fprintf('Starting Newmark for true (damaged) time-history...\n');
for step = 1:nSteps
    tn1 = tvec(step+1);
    % % compute global Kf with theta_true(:,step+1)
    % Kg_time = zeros(Ndof,Ndof);
    % for e = 1:Nelem
    %     dofs = ElementDofs{e};
    %     Ke = ElementKe{e} * theta_true(e, step+1); % scaled elemental stiffness
    %     Kg_time(dofs,dofs) = Kg_time(dofs,dofs) + Ke;
    % end
    % Kf_time = Kg_time(freeDOF, freeDOF);
    % assemble effective stiffness
    Keff = Kf_base*(1-theta_true(step)) + a0*Mf + a1*Cf;
    Fext_n1 = Ff; % static in this example times load factor if needed
    Meff_rhs = Fext_n1 + Mf*(a0*u_true + a2*v_true + a3*a_true) + Cf*(a1*u_true + a4*v_true + a5*a_true);
    % solve
    % guard: Keff may be symmetric pos def -- try chol else backslash
    useChol = false;
    try
        Rchol = chol(Keff,'lower'); useChol = true;
    catch
        useChol = false;
    end
    if useChol
        y = Rchol\Meff_rhs; u_new = Rchol'\y;
    else
        u_new = Keff \ Meff_rhs;
    end
    a_new = a0*(u_new - u_true) - a2*v_true - a3*a_true;
    v_new = v_true + dt*((1-gamma)*a_true + gamma*a_new);
    % store
    Ufree_hist_true(:, step+1) = u_new;
    Vfree_hist_true(:, step+1) = v_new;
    Afree_hist_true(:, step+1) = a_new;
    Utemp = zeros(5*Nnode,1); Utemp(freeDOF) = u_new; Utotal_hist_true(:,step+1) = Utemp;
    % shift
    u_true = u_new; v_true = v_new; a_true = a_new;
    if mod(step, round(nSteps/10))==0, fprintf('  Newmark step %d/%d\n', step, nSteps); end
end
fprintf('True time-history generated.\n');

%% ------------------- Sensor selection --------------------------------
[~, centerNode] = min((mid_coords(:,1)-1).^2 + (mid_coords(:,2)-1).^2);
sensorNodes = centerNode;
% pick two others in interior
candidateNodes = [];
for nn = 1:Nnode
    gDOF = 5*(nn-1)+3;
    if any(freeDOF==gDOF) && nn~=centerNode
        xpos = mid_coords(nn,1); ypos = mid_coords(nn,2);
        if xpos>0.3 && xpos<1.7 && ypos>0.3 && ypos<1.7
            candidateNodes = [candidateNodes; nn];
        end
    end
end
if length(candidateNodes)>=2
    dists = vecnorm(mid_coords(candidateNodes,:) - mid_coords(centerNode,:), 2, 2);
    [~, idx] = max(dists);
    sensorNodes = [sensorNodes; candidateNodes(idx)];
    candidateNodes(idx)=[];
    if ~isempty(candidateNodes), sensorNodes = [sensorNodes; candidateNodes(1)]; end
end
sensorGlobalDOF = 5*(sensorNodes-1)+3;
[~, sensorFreeIdx] = ismember(sensorGlobalDOF, freeDOF);
sensorFreeIdx = sensorFreeIdx(sensorFreeIdx>0);
m = length(sensorFreeIdx);
PHI_s = PHI_red(sensorFreeIdx,:);

fprintf('Sensors chosen at nodes: '); fprintf('%d ', sensorNodes); fprintf('\n');

%% ------------------ Filter parameters & initialization ----------------
% Measurement covariances (regularized)
sigma_meas = 1e-1;
Rk =  1e-2*eye(m);
Rk_UV = 1e-3 * eye(2*m) + 1e-8*eye(2*m);

% state/process noise
Qd_state_mrkf = blkdiag(1e-12*eye(n_modes), 1e-12*eye(n_modes));
Qd_state_dkf  = blkdiag(1e-12*eye(n_modes), 1e-12*eye(n_modes));
% parameter process noise
n_params = 1;
Rk_param = 1e-2*eye(n_params);
param_process_noise = 1e-4;
Qd_param_dkf = param_process_noise * eye(n_params);

% discretize continuous-time modal system correctly (for MRKF)
Z = zeros(n_modes); I_n = eye(n_modes);
A_c = [Z, I_n; -(M_modal\K_modal), -(M_modal\C_modal)];
B_c = [zeros(n_modes); M_modal\eye(n_modes)];
% build Mbig correctly: (2n + n) x (2n + n)
Mbig = [A_c, B_c; zeros(n_modes, 2*n_modes), zeros(n_modes, n_modes)] * dt;
Ebig = expm(Mbig);
Ad = Ebig(1:2*n_modes, 1:2*n_modes);
Bd = Ebig(1:2*n_modes, 2*n_modes+1:end);

% MRKF init
q0 = PHI_red' * Ufree_hist_true(:,1);
qd0 = PHI_red' * Vfree_hist_true(:,1);
x_mrkf = [q0; qd0];
P_mrkf = 1e-4 * eye(2*n_modes);

% DKF init
x_state_dkf = [q0; qd0];
P_state_dkf = 1e-4 * eye(2*n_modes);
x_param_dkf = 0.0;  % initial theta guess (multiplicative on stiffness)
P_param_dkf = 1e-1 * eye(n_params);

% Storage
Qhat_hist_mrkf = zeros(n_modes, nSteps+1);
Qdhat_hist_mrkf = zeros(n_modes, nSteps+1);
Uest_hist_mrkf = zeros(nFree, nSteps+1);
y_meas_hist_mrkf = zeros(m, nSteps+1);
y_hat_hist_mrkf = zeros(m, nSteps+1);

Qhat_hist_dkf = zeros(n_modes, nSteps+1);
Qdhat_hist_dkf = zeros(n_modes, nSteps+1);
Param_hist_dkf = zeros(n_params, nSteps+1);
Uest_hist_dkf = zeros(nFree, nSteps+1);
y_meas_hist_dkf = zeros(m, nSteps+1);
y_hat_hist_dkf = zeros(m, nSteps+1);

Qhat_hist_mrkf(:,1) = q0; Qdhat_hist_mrkf(:,1) = qd0;
Qhat_hist_dkf(:,1) = q0; Qdhat_hist_dkf(:,1) = qd0;
Param_hist_dkf(:,1) = x_param_dkf;
Uest_hist_mrkf(:,1) = PHI_red * q0;
Uest_hist_dkf(:,1) = PHI_red * q0;

U_meas_prev = zeros(m,1); V_meas_prev = zeros(m,1);

% H_linear for MRKF (maps modal states to measured accelerations)
H_linear = PHI_s*[-(M_modal)\K_modal , -(M_modal)\C_modal];

%% ------------------- MRKF & DKF loops --------------------------------
fprintf('Starting MRKF & DKF loops (%d steps)...\n', nSteps);
tic;
eps_reg = 1e-9; % tiny regularizer for matrix inverses
for step = 1:nSteps
    lf = 1; % load factor - here static/1, change if time-dependent forces
    % true measurement at sensors (accelerations) + noise
    a_true_step = Afree_hist_true(:, step+1);
    y_meas_acc = a_true_step(sensorFreeIdx) + sigma_meas * randn(m,1);
    y_meas_hist_mrkf(:, step+1) = y_meas_acc;
    y_meas_hist_dkf(:, step+1) = y_meas_acc;
    % integrate measured acceleration to get pseudo-velocity & displacement measurements
    V_meas_new = V_meas_prev + dt * y_meas_acc;
    U_meas_new = U_meas_prev + dt * V_meas_prev + 0.5 * dt^2 * y_meas_acc;

    %% ---------------- MRKF (baseline) --------------------------------
    u_modal = F_modal * lf;
    x_mrkf_pred = Ad * x_mrkf + Bd * u_modal;
    P_mrkf_pred = Ad * P_mrkf * Ad' + Qd_state_mrkf;
    q_pred = x_mrkf_pred(1:n_modes); qd_pred = x_mrkf_pred(n_modes+1:end);
    qdd_pred = M_modal \ (u_modal - C_modal*qd_pred - K_modal*q_pred);
    y_hat = PHI_s * qdd_pred;
    y_hat_hist_mrkf(:, step+1) = y_hat;
    % MRKF update
    S_mrkf = H_linear * P_mrkf_pred * H_linear' + Rk + eps_reg * eye(size(Rk));
    % robust inverse
    Kk_mrkf = (P_mrkf_pred * H_linear') * pinv(S_mrkf);
    x_mrkf = x_mrkf_pred + Kk_mrkf * (y_meas_acc - y_hat);
    P_mrkf = (eye(2*n_modes) - Kk_mrkf * H_linear) * P_mrkf_pred;
    % store
    Qhat_hist_mrkf(:, step+1) = x_mrkf(1:n_modes);
    Qdhat_hist_mrkf(:, step+1) = x_mrkf(n_modes+1:end);
    Uest_hist_mrkf(:, step+1) = PHI_red * x_mrkf(1:n_modes);

    %% ---------------- DKF (iterative two-filter) ---------------------
    
   
    % Build continuous A_c_ using parameter in modal K = x_param_pred * K_modal
    A_c_ = [zeros(n_modes), eye(n_modes);
            -(M_modal \ ((1-x_param_dkf) * K_modal)), -(M_modal \ C_modal)];
    B_c_ = [zeros(n_modes); M_modal \ eye(n_modes)];
    % Discretize correctly (same pattern)
    Mbig_dkf = [A_c_, B_c_; zeros(n_modes, 2*n_modes), zeros(n_modes, n_modes)] * dt;
    Ebig_dkf = expm(Mbig_dkf);
    Ad_dkf = Ebig_dkf(1:2*n_modes, 1:2*n_modes);
    Bd_dkf = Ebig_dkf(1:2*n_modes, 2*n_modes+1:end);

    % State prediction (DKF)
    u_modal_est = F_modal * lf; % assume known input (no param dependence here)
    x_state_pred = Ad_dkf * x_state_dkf + Bd_dkf * u_modal_est;
    P_state_pred = Ad_dkf * P_state_dkf * Ad_dkf' + Qd_state_dkf;

    % Predicted measurement (accelerations) under predicted state & param
    q_pred = x_state_pred(1:n_modes); qd_pred = x_state_pred(n_modes+1:end);
    qdd_pred = M_modal \ (u_modal_est - C_modal * qd_pred - ((1-x_param_dkf)* K_modal) * q_pred);
    y_hat_dkf = PHI_s * qdd_pred;
    y_hat_hist_dkf(:, step+1) = y_hat_dkf;

   

    % Build H_linear_dkf mapping x_state to [U;V] measured (2m)
    % we measure U and V derived from integration of accel; the relation used in the script:
    % U = PHI_s * q  (m x n_modes) times q -> measured displacement at sensors
    % V = PHI_s * qd
    % H_U = PHI_s; H_V = PHI_s;
    % H_linear_dkf = [H_U, zeros(m, n_modes); zeros(m,n_modes), H_V]; % (2m x 2n_modes)
    H_linear_dkf = PHI_s*[ -M_modal\((1-x_param_dkf)*K_modal),   -(M_modal\C_modal)];
    D_state_dkf = PHI_s*(M_modal\u_modal_est);
    % Now measurement residual for [U;V]
    y_state_residual = y_meas_acc - (H_linear_dkf*x_state_pred + D_state_dkf);

    % Compute S_state robustly
    S_state = H_linear_dkf * P_state_pred * H_linear_dkf' + Rk ;
    % safeguard: if S_state not full rank, use pinv
    Kk_state = (P_state_pred * H_linear_dkf') * pinv(S_state);
    % update the state
    x_state_dkf = x_state_pred + Kk_state * y_state_residual;
    P_state_dkf = (eye(2*n_modes) - Kk_state * H_linear_dkf) * P_state_pred;

%% ------- Parameter ------------

    % Parameter prediction (random walk)
    x_param_pred = x_param_dkf + 1e-12;
    P_param_pred = P_param_dkf + Qd_param_dkf;

     % Measurement sensitivity wrt theta:
    H_param =  PHI_s * (M_modal \ (K_modal * x_state_dkf(1:n_modes))); % m x n_params (here n_params=1)
    D_param = PHI_s * (M_modal \(u_modal_est - K_modal*x_state_dkf(1:n_modes) - C_modal*x_state_dkf(n_modes+1:end)));
    % Parameter update using acceleration residual (scalar/m vector)
    y_residual_param = y_meas_acc - (H_param*x_param_pred + D_param);
    S_param = H_param * P_param_pred * H_param' + Rk_param ;
    Kk_param = (P_param_pred * H_param') * pinv(S_param);
    x_param_dkf = x_param_pred + Kk_param * y_residual_param;
    P_param_dkf = (eye(n_params) - Kk_param * H_param) * P_param_pred;

    % store
    Qhat_hist_dkf(:, step+1) = x_state_dkf(1:n_modes);
    Qdhat_hist_dkf(:, step+1) = x_state_dkf(n_modes+1:end);
    Param_hist_dkf(:, step+1) = x_param_dkf;
    Uest_hist_dkf(:, step+1) = PHI_red * x_state_dkf(1:n_modes);
    y_meas_hist_dkf(:, step+1) = y_meas_acc;

    % update measured-integration states
    U_meas_prev = U_meas_new; V_meas_prev = V_meas_new;

    % diagnostics printing
    if mod(step, max(1,round(nSteps/10)))==0
        fprintf('Step %d/%d: theta_est=%.4f\n', step, nSteps, x_param_dkf);
    end

end
toc;

%% ------------------- Postprocessing & plots ---------------------------
% Convert estimates to full DOF
Utotal_hist_mrkf = zeros(5*Nnode, nSteps+1);
Utotal_hist_dkf = zeros(5*Nnode, nSteps+1);
for i = 1:nSteps+1
    temp = zeros(5*Nnode,1); temp(freeDOF) = Uest_hist_mrkf(:,i); Utotal_hist_mrkf(:,i) = temp;
    temp2 = zeros(5*Nnode,1); temp2(freeDOF) = Uest_hist_dkf(:,i); Utotal_hist_dkf(:,i) = temp2;
end

% Plot comparison at sensor nodes
figure('Position',[100 100 1200 800]);
for s = 1:m
    subplot(m,1,s);
    sensor_node = sensorNodes(s);
    sensor_dof = 5*(sensor_node-1)+3;
    true_w = squeeze(Utotal_hist_true(sensor_dof,:));
    mrkf_w = squeeze(Utotal_hist_mrkf(sensor_dof,:));
    dkf_w = squeeze(Utotal_hist_dkf(sensor_dof,:));
    plot(tvec, true_w,'k-','LineWidth',1.5); hold on;
    plot(tvec, mrkf_w,'g--'); plot(tvec, dkf_w,'b-.'); grid on;
    title(sprintf('Sensor %d (Node %d) displacement w', s, sensor_node));
    legend('True','MRKF','DKF');
end

% Parameter plot
figure;
plot(tvec, Param_hist_dkf(1,:),'r-','LineWidth',1.6); hold on;
% plot true mean theta for damaged region (not exactly scalar if multiple elements)
%theta_mean_true = mean(theta_true(damaged_elems,:),1); % average of damaged elements
plot(tvec, theta_true,'k--','LineWidth',1.2);
xlabel('Time (s)'); ylabel('Estimated \theta (DKF)'); legend('DKF estimate','True avg \theta');

fprintf('Finished. Final estimated theta = %.6f\n', Param_hist_dkf(1,end));
