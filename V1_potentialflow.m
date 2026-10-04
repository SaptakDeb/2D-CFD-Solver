clear; clc; close all;

%% Grid
Nx = 160;
Ny = 100;

Lx = 4;
Ly = 2;

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);

x = linspace(-1,Lx-1,Nx);
y = linspace(-Ly/2,Ly/2,Ny);

[X,Y] = meshgrid(x,y);

%% Flow conditions
gamma = 1.4;
R = 287;

M_inf = 0.3;
T_inf = 300;
p_inf = 101325;

rho_inf = p_inf/(R*T_inf);
a_inf = sqrt(gamma*R*T_inf);
U_inf = M_inf*a_inf;

cp = gamma*R/(gamma-1);

% Stagnation speed of sound and speed limit (keeps the density law valid
% when the sharp corners produce very large discrete velocities)
a0_sq = a_inf^2 + (gamma-1)/2*U_inf^2;
V2_max = 0.9*2*a0_sq/(gamma+1);

%% Rectangular solid
x1 = 0.8;
x2 = 1.8;

y1 = -0.35;
y2 = 0.35;

solid = (X >= x1 & X <= x2 & Y >= y1 & Y <= y2);
fluid = ~solid;

%% Potential flow initial condition
phi = U_inf*X;

% Density ratio rho/rho_inf (updated every iteration)
rho_r = ones(Ny,Nx);

%% Solver settings
omega = 1.85;        % SOR factor for phi
omega_rho = 0.8;     % under-relaxation for density
tol = 1e-6;          % scaled by U_inf*dx
maxIter = 5000;
printEvery = 100;

%% Index sets (interior nodes, red-black ordering)
J = 2:Ny-1;
I = 2:Nx-1;

[II,JJ] = meshgrid(I,J);

red = (mod(II+JJ,2) == 0);
black = ~red;

fluid_in = fluid(J,I);

% Neighbour solid flags -> zero mass flux through that face
fE = double(~solid(J,I+1));
fW = double(~solid(J,I-1));
fN = double(~solid(J+1,I));
fS = double(~solid(J-1,I));

%% Solver start
fprintf('\n');
fprintf('Full Potential Solver\n');
fprintf('Equation          : d(rho*phi_x)/dx + d(rho*phi_y)/dy = 0\n');
fprintf('Scheme            : conservative finite volume, red-black SOR\n');
fprintf('Grid              : %d x %d\n',Nx,Ny);
fprintf('Mach number       : %.3f\n',M_inf);
fprintf('SOR factor        : %.2f\n',omega);
fprintf('Density relaxation: %.2f\n',omega_rho);
fprintf('Tolerance         : %.1e\n',tol);
fprintf('Maximum iterations: %d\n',maxIter);
fprintf('\n');

errHist = nan(maxIter,1);

tic;

%% Red-black SOR iteration
for iter = 1:maxIter

    phi_old = phi;

    % Face densities (arithmetic average of neighbouring nodes)
    rE = 0.5*(rho_r(J,I) + rho_r(J,I+1));
    rW = 0.5*(rho_r(J,I) + rho_r(J,I-1));
    rN = 0.5*(rho_r(J,I) + rho_r(J+1,I));
    rS = 0.5*(rho_r(J,I) + rho_r(J-1,I));

    % Face coefficients, solid faces carry no flux (dphi/dn = 0)
    aE = rE.*fE/dx^2;
    aW = rW.*fW/dx^2;
    aN = rN.*fN/dy^2;
    aS = rS.*fS/dy^2;

    aP = aE + aW + aN + aS;
    aP(aP == 0) = 1;

    % Two colour sweeps, black nodes use the updated red values
    for color = 1:2

        if color == 1
            mask = red & fluid_in;
        else
            mask = black & fluid_in;
        end

        numerator = aE.*phi(J,I+1) + aW.*phi(J,I-1) ...
                  + aN.*phi(J+1,I) + aS.*phi(J-1,I);

        phi_gs = numerator./aP;

        phi(J,I) = phi(J,I) + omega*(phi_gs - phi(J,I)).*mask;

        % Top/bottom: zero normal velocity
        phi(1,:)  = phi(2,:);
        phi(Ny,:) = phi(Ny-1,:);

    end

    % Inlet/outlet: phi = U_inf*x (never updated, kept fixed)
    phi(:,1)  = U_inf*X(:,1);
    phi(:,Nx) = U_inf*X(:,Nx);

    % Error (scaled so it is independent of grid size and speed)
    err = max(max(abs(phi - phi_old)))/(U_inf*dx);
    errHist(iter) = err;

    % Divergence check
    if ~isfinite(err)
        fprintf('\n');
        fprintf('SOLVER DIVERGED\n');
        fprintf('Iteration : %d\n',iter);
        fprintf('Error     : %g\n',err);
        fprintf('Time      : %.2f s\n',toc);
        break;
    end

    % Velocity at nodes (ghost value = current value next to a solid face)
    [u,v] = gradient(phi,dx,dy);

    phiC = phi(J,I);
    phiE = phi(J,I+1);
    phiW = phi(J,I-1);
    phiN = phi(J+1,I);
    phiS = phi(J-1,I);

    phiE(~logical(fE)) = phiC(~logical(fE));
    phiW(~logical(fW)) = phiC(~logical(fW));
    phiN(~logical(fN)) = phiC(~logical(fN));
    phiS(~logical(fS)) = phiC(~logical(fS));

    u(J,I) = (phiE - phiW)/(2*dx);
    v(J,I) = (phiN - phiS)/(2*dy);

    V2 = u.^2 + v.^2;
    V2c = min(V2,V2_max);

    % Isentropic density update
    T_r = 1 + (gamma-1)/2*M_inf^2*(1 - V2c/U_inf^2);
    rho_new = T_r.^(1/(gamma-1));

    rho_r = (1-omega_rho)*rho_r + omega_rho*rho_new;

    % Live progress
    if mod(iter,printEvery) == 0
        elapsed = toc;
        percent = 100*iter/maxIter;

        fprintf(['Iteration: %5d / %5d   ' ...
                 'Progress: %6.2f %%   ' ...
                 'Error: %.3e   ' ...
                 'Time: %.1f s\n'], ...
                 iter,maxIter,percent,err,elapsed);
    end

    % Convergence check
    if err < tol
        elapsed = toc;
        fprintf('\n');
        fprintf('Converged\n');
        fprintf('Iterations : %d\n',iter);
        fprintf('Final error: %.3e\n',err);
        fprintf('Time       : %.2f s\n',elapsed);
        break;
    end

end

errHist = errHist(1:iter);

%% Maximum iteration warning
if iter == maxIter && err >= tol
    fprintf('\n');
    fprintf('Maximum iterations reached\n');
    fprintf('Final error: %.3e\n',err);
    fprintf('Time       : %.2f s\n',toc);
end

%% Velocity
fprintf('\nCalculating velocity...\n');

V2 = u.^2 + v.^2;
V = sqrt(V2);

nClipped = sum(V2(fluid) > V2_max);

%% Conservation check (mass flux through inlet and outlet)
mass_in  = trapz(y(:), rho_r(:,2).*u(:,2))*rho_inf;
mass_out = trapz(y(:), rho_r(:,Nx-1).*u(:,Nx-1))*rho_inf;

% Mask solid
u(solid) = NaN;
v(solid) = NaN;
V(solid) = NaN;
V2(solid) = NaN;

%% Compressible temperature
fprintf('Calculating temperature...\n');

T = T_inf + (U_inf^2 - V2)/(2*cp);

% Prevent nonphysical negative/zero temperatures
T(T <= 0) = NaN;
T(solid) = NaN;

%% Density
rho = rho_inf .* (T/T_inf).^(1/(gamma-1));
rho(solid) = NaN;

%% Pressure
fprintf('Calculating pressure...\n');

p = rho .* R .* T;
p(solid) = NaN;

%% Mach number
fprintf('Calculating Mach number...\n');

a = sqrt(gamma*R*T);
Mach = V./a;
Mach(solid) = NaN;

%% Final values
fprintf('\n');
fprintf('Solution complete\n');
fprintf('Freestream velocity : %.3f m/s\n',U_inf);
fprintf('Maximum velocity    : %.3f m/s\n',max(V(:)));
fprintf('Minimum pressure    : %.2f Pa\n',min(p(:)));
fprintf('Maximum pressure    : %.2f Pa\n',max(p(:)));
fprintf('Maximum Mach number : %.4f\n',max(Mach(~isnan(Mach))));
fprintf('Clipped nodes       : %d\n',nClipped);
fprintf('Mass flux in        : %.4f kg/s per m\n',mass_in);
fprintf('Mass flux out       : %.4f kg/s per m\n',mass_out);
fprintf('Mass flux imbalance : %.3f %%\n',100*(mass_out-mass_in)/mass_in);

%% Plot 1: velocity magnitude
figure;
contourf(X,Y,V,30,'LineColor','none');
hold on;
contour(X,Y,double(solid),[1 1],'k','LineWidth',2);
axis equal;
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
colorbar;
xlabel('x'); ylabel('y');
title('Velocity Magnitude');

%% Plot 2: pressure
figure;
contourf(X,Y,p,30,'LineColor','none');
hold on;
contour(X,Y,double(solid),[1 1],'k','LineWidth',2);
axis equal;
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
colorbar;
xlabel('x'); ylabel('y');
title('Pressure');

%% Plot 3: Mach number
figure;
contourf(X,Y,Mach,30,'LineColor','none');
hold on;
contour(X,Y,double(solid),[1 1],'k','LineWidth',2);
axis equal;
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
colorbar;
xlabel('x'); ylabel('y');
title('Mach Number');

%% Plot 4: velocity field
figure;
skip = 4;
quiver( ...
    X(1:skip:end,1:skip:end), ...
    Y(1:skip:end,1:skip:end), ...
    u(1:skip:end,1:skip:end), ...
    v(1:skip:end,1:skip:end));
hold on;
contour(X,Y,double(solid),[1 1],'k','LineWidth',2);
axis equal;
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
xlabel('x'); ylabel('y');
title('Velocity Field');

%% Plot 5: convergence history
figure;
semilogy(1:iter,errHist,'LineWidth',1.5);
grid on;
xlabel('Iteration'); ylabel('Scaled error');
title('Convergence History');

fprintf('\nAll plots generated.\n');
