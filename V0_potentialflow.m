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

%% Rectangular solid
x1 = 0.8;
x2 = 1.8;

y1 = -0.35;
y2 = 0.35;

solid = (X >= x1 & X <= x2 & Y >= y1 & Y <= y2);

%% Potential flow initial condition
phi = U_inf*X;

%% Solver settings
omega = 0.8;
tol = 1e-7;
maxIter = 100;

%% Solver start
fprintf('\n');
fprintf('Potential Flow Solver\n');
fprintf('Grid              : %d x %d\n',Nx,Ny);
fprintf('Mach number       : %.3f\n',M_inf);
fprintf('Relaxation factor : %.2f\n',omega);
fprintf('Tolerance         : %.1e\n',tol);
fprintf('Maximum iterations: %d\n',maxIter);
fprintf('\n');

tic;

%% Gauss-Seidel iteration
for iter = 1:maxIter

    phi_old = phi;

    % Sweep interior fluid nodes
    for j = 2:Ny-1
        for i = 2:Nx-1

            % Skip solid cells
            if solid(j,i)
                continue;
            end

            % Finite-difference weights
            wE = 1/dx^2;
            wW = 1/dx^2;
            wN = 1/dy^2;
            wS = 1/dy^2;

            numerator = 0;
            denominator = 0;

            % East
            if solid(j,i+1)
                % zero normal velocity -> ghost value = current value
                denominator = denominator + wE;
            else
                numerator = numerator + wE*phi(j,i+1);
                denominator = denominator + wE;
            end

            % West
            if solid(j,i-1)
                denominator = denominator + wW;
            else
                numerator = numerator + wW*phi(j,i-1);
                denominator = denominator + wW;
            end

            % North
            if solid(j+1,i)
                denominator = denominator + wN;
            else
                numerator = numerator + wN*phi(j+1,i);
                denominator = denominator + wN;
            end

            % South
            if solid(j-1,i)
                denominator = denominator + wS;
            else
                numerator = numerator + wS*phi(j-1,i);
                denominator = denominator + wS;
            end

            % Gauss-Seidel update
            phi_new = numerator/denominator;
            phi(j,i) = (1-omega)*phi(j,i) + omega*phi_new;

        end
    end

    % Far-field boundaries: phi = U_inf*x
    phi(:,1)  = U_inf*X(:,1);
    phi(:,Nx) = U_inf*X(:,Nx);
    phi(1,:)  = U_inf*X(1,:);
    phi(Ny,:) = U_inf*X(Ny,:);

    % Error
    err = max(abs(phi(:)-phi_old(:)));

    % Divergence check
    if ~isfinite(err)
        fprintf('\n');
        fprintf('SOLVER DIVERGED\n');
        fprintf('Iteration : %d\n',iter);
        fprintf('Error     : %g\n',err);
        fprintf('Time      : %.2f s\n',toc);
        break;
    end

    % Live progress
    if mod(iter,10) == 0
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

%% Maximum iteration warning
if iter == maxIter && err >= tol
    fprintf('\n');
    fprintf('Maximum iterations reached\n');
    fprintf('Final error: %.3e\n',err);
    fprintf('Time       : %.2f s\n',toc);
end

%% Velocity
fprintf('\nCalculating velocity...\n');

[u,v] = gradient(phi,dx,dy);

V = sqrt(u.^2 + v.^2);
V2 = u.^2 + v.^2;

% Mask solid
u(solid) = NaN;
v(solid) = NaN;
V(solid) = NaN;
V2(solid) = NaN;

%% Compressible temperature
fprintf('Calculating temperature...\n');

cp = gamma*R/(gamma-1);
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
fprintf('Maximum Mach number : %.4f\n',max(Mach(:),[],'omitnan'));

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

fprintf('\nAll plots generated.\n');
