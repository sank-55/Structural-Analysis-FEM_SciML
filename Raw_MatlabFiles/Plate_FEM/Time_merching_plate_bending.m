%% Plate_Kirchhoff_FEM_NewmarkPlots.m
% Full Kirchhoff plate FEM with Newmark-beta time marching and
% the animation + center-node time-history plots you requested.
% 3 DOF per node (w, theta_x, theta_y)

clear; close all; clc;

%% ---------------- MATERIAL & GEOMETRY ----------------
E = 210e9;
h = 0.005;
rho = 8050;
nu = 0.3;

Lx = 2; Ly = 2;
q0 = -1500;        % base distributed transverse load (N/m^2)

% mesh (dx = dy)
nx = 20; ny = 20;
dx = Lx/nx; dy = Ly/ny;
Nelem = nx*ny;
Nnode = (nx+1)*(ny+1);
Ndof = 3*Nnode;    % [w, theta_x, theta_y] per node

%% ---------------- BENDING D MATRIX ----------------
D = (E*h^3 / (12*(1 - nu^2))) * [1, nu, 0;
                                 nu, 1, 0;
                                 0, 0, (1 - nu)/2];

%% ---------------- NODES & CONNECTIVITY ----------------
coords = zeros(Nnode,2);
ncount = 1;
for j = 0:ny
    for i = 0:nx
        coords(ncount,:) = [i*dx, j*dy];
        ncount = ncount + 1;
    end
end

c = zeros(Nelem,4);
elem = 0;
for j = 1:ny
    for i = 1:nx
        n1 = (j-1)*(nx+1) + i;
        n2 = n1 + 1;
        n3 = n2 + (nx+1);
        n4 = n1 + (nx+1);
        elem = elem + 1;
        c(elem,:) = [n1 n2 n3 n4];
    end
end

%% ---------------- GAUSS POINTS & SHAPE FUNCTIONS (symbolic) ----------------
gp = [-1/sqrt(3), 1/sqrt(3)];
wgp = [1,1];

syms xi eta real
le = 1;
H1 = 1 - 3*xi^2 + 2*xi^3;
H2 = le*(xi - 2*xi^2 + xi^3);
H3 = 3*xi^2 - 2*xi^3;
H4 = le*(-xi^2 + xi^3);
H_xi = [H1; H2; H3; H4];
H_eta = subs(H_xi, xi, eta);

% 12 shape functions (same order you used)
N1  = H_xi(1) * H_eta(1);
N2  = H_xi(2) * H_eta(1);
N3  = H_xi(1) * H_eta(2);
N4  = H_xi(3) * H_eta(1);
N5  = H_xi(4) * H_eta(1);
N6  = H_xi(3) * H_eta(2);
N7  = H_xi(3) * H_eta(3);
N8  = H_xi(4) * H_eta(3);
N9  = H_xi(3) * H_eta(4);
N10 = H_xi(1) * H_eta(3);
N11 = H_xi(2) * H_eta(3);
N12 = H_xi(1) * H_eta(4);

N_sym = [N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12];

N_xixi  = simplify(diff(N_sym, xi, 2));
N_etaeta = simplify(diff(N_sym, eta, 2));
N_xieta = simplify(diff(diff(N_sym, xi), eta));

% For linear mapping xi->x, eta->y: J = [2/dx 0; 0 2/dy] -> detJ = 4/(dx*dy)
detJ_const = (2/dx)*(2/dy);

% Precompute at 4 Gauss points
numGP = 4;
N_at_gp = zeros(numGP, 12);
B_at_gp = zeros(3,12,numGP);
wts_at_gp = zeros(numGP,1);
cnt = 0;
for i = 1:2
    for j = 1:2
        cnt = cnt + 1;
        xi_val = gp(i); eta_val = gp(j);
        Nval = double(subs(N_sym, [xi,eta], [xi_val, eta_val]));
        N_xixi_val = double(subs(N_xixi, [xi,eta], [xi_val, eta_val]));
        N_etaeta_val = double(subs(N_etaeta, [xi,eta], [xi_val, eta_val]));
        N_xieta_val = double(subs(N_xieta, [xi,eta], [xi_val, eta_val]));
        Bval = (4/dx^2) * [N_xixi_val; N_etaeta_val; 2*N_xieta_val];
        N_at_gp(cnt,:) = Nval;
        B_at_gp(:,:,cnt) = Bval;
        wts_at_gp(cnt) = wgp(i)*wgp(j);
    end
end

%% ---------------- ASSEMBLY ----------------
Kg = sparse(Ndof, Ndof);
Mg = sparse(Ndof, Ndof);
Fg = zeros(Ndof,1);

for e = 1:Nelem
    nodes = c(e,:);
    Ke = zeros(12,12);
    Me = zeros(12,12);
    Fe = zeros(12,1);
    for gpidx = 1:numGP
        Nval = N_at_gp(gpidx,:)';
        Bval = B_at_gp(:,:,gpidx);
        wt = wts_at_gp(gpidx);
        Ke = Ke + (Bval')*D*Bval*detJ_const*wt;
        Me = Me + (rho*h) * (Nval * Nval') * detJ_const * wt;
        Fe = Fe + q0 * detJ_const * wt * Nval;
    end
    dofs = reshape([3*nodes-2; 3*nodes-1; 3*nodes],1,[]);
    Kg(dofs,dofs) = Kg(dofs,dofs) + Ke;
    Mg(dofs,dofs) = Mg(dofs,dofs) + Me;
    Fg(dofs) = Fg(dofs) + Fe;
end

%% ---------------- BOUNDARY CONDITIONS (clamp left edge like your original) ----------------
tol = 1e-9;
left_nodes = find(abs(coords(:,1)) < tol);
fixedDOF = reshape([3*left_nodes-2; 3*left_nodes-1; 3*left_nodes],1,[]);
freeDOF = setdiff(1:Ndof, fixedDOF);

Kf = Kg(freeDOF, freeDOF);
Mf = Mg(freeDOF, freeDOF);
Ff = Fg(freeDOF);

%% ---------------- Modal extraction for Rayleigh damping ----------------
nEig = 10;
try
    [PHI_all, W2_all] = eigs(Kf, Mf, min(nEig, size(Kf,1)), 'smallestabs');
    omega_all = sqrt(abs(diag(W2_all)));
catch
    [PHIs, W2s] = eig(full(Kf), full(Mf));
    omega_all = sqrt(abs(diag(W2s)));
    [omega_all, idxs] = sort(omega_all);
    PHI_all = PHIs(:, idxs(1:min(nEig,length(idxs))));
end

if isempty(omega_all)
    error('Modal extraction failed — check Kf/Mf.');
end

% Rayleigh damping (two-point fit)
zeta_target = 0.02;
w1 = omega_all(1);
w2 = omega_all(min(3,length(omega_all)));
% avoid division by zero
if w1==0 || w2==0
    alphaR = 0; betaR = 0;
else
    A = [1/w1, w1; 1/w2, w2];
    bvec = 2*zeta_target * [1;1];
    ab = A\bvec;
    alphaR = ab(1); betaR = ab(2);
end
Cf = alphaR * Mf + betaR * Kf;
fprintf('Rayleigh: alpha=%.4e, beta=%.4e (w1=%.3f, w2=%.3f)\n', alphaR, betaR, w1, w2);

%% ---------------- Newmark parameters & time stepping setup ----------------
beta_nm = 1/4; gamma_nm = 1/2;  % average-acceleration

dt = 0.001;
Ttotal = 1.0;
tvec = 0:dt:Ttotal;
nSteps = length(tvec)-1;

Kred = Kf; Mred = Mf; Cred = Cf;
F_static = Ff;

% ramped load to avoid shock (you used tramp earlier)
tramp = 0.01;
load_factor = @(tt) min(1, tt/tramp);

% Newmark constants
a0 = 1/(beta_nm*dt^2);
a1 = gamma_nm/(beta_nm*dt);
a2 = 1/(beta_nm*dt);
a3 = 1/(2*beta_nm)-1;
a4 = gamma_nm/beta_nm - 1;
a5 = dt*(gamma_nm/(2*beta_nm) - 1);

Keff = Kred + a0*Mred + a1*Cred;

% try Cholesky factorization for speed; fallback to backslash
useChol = false;
Rchol = [];
try
    Rchol = chol(Keff);
    useChol = true;
catch
    useChol = false;
end

% initial conditions
nFree = size(Kred,1);
u = zeros(nFree,1);
v = zeros(nFree,1);
a = Mred \ ( (F_static * load_factor(0)) - Cred*v - Kred*u );

% Storage
Ufree_hist = zeros(nFree, nSteps+1);
Ufull_hist = zeros(Ndof, nSteps+1);
w_time = zeros(Nnode, nSteps+1);

% initial
Ufree_hist(:,1) = u;
Ufull = zeros(Ndof,1);
Ufull(freeDOF) = u;
Ufull_hist(:,1) = Ufull;
w_time(:,1) = Ufull(1:3:end);

fprintf('Starting Newmark integration: %d steps, dt=%.4g s\n', nSteps, dt);
tic;
for step = 1:nSteps
    tn1 = tvec(step+1);
    Fext_n1 = F_static * load_factor(tn1);
    Meff_rhs = Fext_n1 + Mred*( a0*u + a2*v + a3*a ) + Cred*( a1*u + a4*v + a5*a );
    if useChol
        % solve Keff * u_new = Meff_rhs via Cholesky
        y = Rchol'\Meff_rhs;
        u_new = Rchol\y;
    else
        u_new = Keff \ Meff_rhs;
    end
    a_new = a0*(u_new - u) - a2*v - a3*a;
    v_new = v + dt*((1-gamma_nm)*a + gamma_nm*a_new);
    % store
    Ufree_hist(:, step+1) = u_new;
    Ufull = zeros(Ndof,1);
    Ufull(freeDOF) = u_new;
    Ufull_hist(:, step+1) = Ufull;
    w_time(:, step+1) = Ufull(1:3:end);
    % shift
    u = u_new; v = v_new; a = a_new;
    if mod(step, max(1,round(nSteps/10)))==0
        fprintf('  step %d/%d, t=%.4f s\n', step, nSteps, tn1);
    end
end
toc;
fprintf('Time integration finished.\n');

%% ---------------- Post-processing: center node time history ----------------
% find center node (closest to (Lx/2, Ly/2))
[~, centerIdx] = min( (coords(:,1)-Lx/2).^2 + (coords(:,2)-Ly/2).^2 );
w_center = w_time(centerIdx, :);

figure;
plot(tvec, w_center, 'LineWidth', 1.4);
xlabel('Time (s)'); ylabel('w at center (m)');
title('Transverse displacement at center node vs time');
grid on;

%% ---------------- Animation of mid-surface deflection ----------------
% Amplification parameters
AMP_MULTIPLIER = 2;   % tweak if you want larger visual amplification

x = coords(:,1); y = coords(:,2);
xq = linspace(min(x), max(x), 80);
yq = linspace(min(y), max(y), 80);
[Xq, Yq] = meshgrid(xq, yq);

W_all = w_time;  % Nnode x time
maxAbsW = max(abs(W_all(:)));
if maxAbsW == 0
    warning('All deflections are zero -> check loading or BCs.');
    maxAbsW = 1e-12;
end

geomScale = max(max(xq)-min(xq), max(yq)-min(yq));
targetFrac = 0.05;
base_amp = (geomScale * targetFrac) / maxAbsW;
base_amp = max(base_amp, 1);
amp = base_amp * AMP_MULTIPLIER;
amp = min(amp, 1e8);

fprintf('Animation amp: base=%.3g, multiplier=%g, final amp=%.3g\n', base_amp, AMP_MULTIPLIER, amp);

% prepare interpolant
Finterp = scatteredInterpolant(x, y, W_all(:,1), 'natural', 'none');
Zq0 = Finterp(Xq, Yq);

figure('Units','normalized','Position',[0.05 0.05 0.6 0.7]);
hSurf = surf(Xq, Yq, amp * Zq0);
shading interp; axis equal; view(45,30);
set(hSurf, 'EdgeColor','none');
xlabel('x'); ylabel('y'); zlabel('w_{amp}');
colormap(parula); colorbar;
caxis(amp*[-maxAbsW, maxAbsW]);

% time-history small plot
figure('Units','normalized','Position',[0.7 0.65 0.25 0.25]);
plot(tvec, w_center, '-', 'LineWidth', 1.2); hold on;
hNow = plot(tvec(1), w_center(1), 'or', 'MarkerFaceColor','r');
xlabel('Time (s)'); ylabel('w at center (m)'); title('Center node time history'); grid on;

% animation loop
frameStep = max(1, round(nSteps/100));
for k = 1:frameStep:(nSteps+1)
    Finterp.Values = W_all(:,k);
    Zq_real = Finterp(Xq, Yq);
    set(hSurf, 'ZData', amp * Zq_real);
    title(sprintf('Amplified deflection (x%.3g), t = %.4f s', amp, tvec(k)));
    % update time-history marker
    set(0,'CurrentFigure',findobj('Type','figure','Name','') ); % keep safe (no-op) 
    set( findobj('Type','line','Marker','o'),'XData', tvec(k), 'YData', w_center(k)); %#ok<SAGROW>
    drawnow;
end

%% ---------------- Save or return results if desired ----------------
% Example: save center time history
% save('center_deflection.mat', 'tvec', 'w_center');

% END of script
