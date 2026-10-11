clear; clc; close all;

%% Flow conditions
gamma = 1.4;
R = 287;

M_inf = 0.3;
T_inf = 300;
p_inf = 101325;

aoa_deg = 0;                 % angle of attack [deg], positive = nose up
aoa = aoa_deg*pi/180;

rho_inf = p_inf/(R*T_inf);
a_inf = sqrt(gamma*R*T_inf);
U_inf = M_inf*a_inf;
q_inf = 0.5*rho_inf*U_inf^2;
cp = gamma*R/(gamma-1);

%% Rectangular solid (rotated about its own centre)
Lb = 1.0;                    % length (chord) of the rectangle
Hb = 0.7;                    % thickness of the rectangle
xc_b = 1.3;                  % rectangle centre x (same as v1: 0.8 to 1.8)
yc_b = 0.0;                  % rectangle centre y

%% Domain size
x_min = xc_b - 15;           % upstream extent
x_max = xc_b + 25;           % downstream extent (longer, for the wake)
y_min = -15;
y_max =  15;

%% Grid (uniform refined box around the body, geometric stretching outside)
h0    = 0.025;               % cell size in the refined box
r     = 1.06;                % stretching ratio outside the box
a_ref = 1.0;                 % half-width of the refined box around the body

% x direction
nRef = round(2*a_ref/h0);
xin  = linspace(xc_b-a_ref, xc_b+a_ref, nRef+1);
LL   = (xc_b-a_ref) - x_min;
LR   = x_max - (xc_b+a_ref);
nL   = ceil(log(1 + LL*(r-1)/h0)/log(r));
nR   = ceil(log(1 + LR*(r-1)/h0)/log(r));
wL   = h0*r.^(1:nL);  wL = wL*(LL/sum(wL));
wR   = h0*r.^(1:nR);  wR = wR*(LR/sum(wR));
xf   = [fliplr((xc_b-a_ref) - cumsum(wL)), xin, (xc_b+a_ref) + cumsum(wR)];

% y direction
yin  = linspace(yc_b-a_ref, yc_b+a_ref, nRef+1);
LB   = (yc_b-a_ref) - y_min;
LT   = y_max - (yc_b+a_ref);
nB   = ceil(log(1 + LB*(r-1)/h0)/log(r));
nT   = ceil(log(1 + LT*(r-1)/h0)/log(r));
wB   = h0*r.^(1:nB);  wB = wB*(LB/sum(wB));
wT   = h0*r.^(1:nT);  wT = wT*(LT/sum(wT));
yf   = [fliplr((yc_b-a_ref) - cumsum(wB)), yin, (yc_b+a_ref) + cumsum(wT)];

x  = 0.5*(xf(1:end-1) + xf(2:end));
y  = 0.5*(yf(1:end-1) + yf(2:end));
dx = diff(xf);
dy = diff(yf);

Nx = numel(x);
Ny = numel(y);

[X,Y] = meshgrid(x,y);
DX = ones(Ny,1)*dx;
DY = dy(:)*ones(1,Nx);

% Padded (ghost cell) widths and centres used by the reconstruction
dxp = [dx(1), dx, dx(end)];
dyp = [dy(1), dy, dy(end)];
xcp = [x(1)-dx(1), x, x(end)+dx(end)];
ycp = [y(1)-dy(1), y, y(end)+dy(end)];

%% Rotate the geometry (clockwise by aoa, so positive aoa = nose up)
ca = cos(aoa);
sa = sin(aoa);

xb =  (X-xc_b)*ca - (Y-yc_b)*sa;     % cell centres in body-fixed frame
yb =  (X-xc_b)*sa + (Y-yc_b)*ca;

solid = (abs(xb) <= Lb/2 & abs(yb) <= Hb/2);
fluid = ~solid;

Sp = false(Ny+2,Nx+2);
Sp(2:Ny+1,2:Nx+1) = solid;

% Weights for the pressure force on the (staircase) body surface
sE = double(Sp(2:Ny+1,3:Nx+2));
sW = double(Sp(2:Ny+1,1:Nx));
sN = double(Sp(3:Ny+2,2:Nx+1));
sS = double(Sp(1:Ny,2:Nx+1));
wFx = double(fluid).*(sE - sW).*DY;
wFy = double(fluid).*(sN - sS).*DX;

%% Solver settings
order     = 2;               % 1 = first order, 2 = MUSCL (minmod)
nFirst    = 300;             % first-order start-up iterations
CFL       = 0.8;
rampIter  = 300;             % CFL ramp length
maxIter   = 5000;
tol       = 1e-4;            % density residual drop (relative to first iteration)
printEvery = 2;

rho_floor = 0.01*rho_inf;
p_floor   = 0.01*p_inf;

%% Initial condition: uniform freestream
Uc = zeros(Ny,Nx,4);
Uc(:,:,1) = rho_inf;
Uc(:,:,2) = rho_inf*U_inf;
Uc(:,:,3) = 0;
Uc(:,:,4) = p_inf/(gamma-1) + 0.5*rho_inf*U_inf^2;

Wp = zeros(Ny+2,Nx+2,4);     % padded primitive variables [rho u v p]
Wp(:,:,1) = rho_inf;
Wp(:,:,2) = U_inf;
Wp(:,:,3) = 0;
Wp(:,:,4) = p_inf;

%% Solver start
fprintf('\n');
fprintf('2D Euler Solver\n');
fprintf('Equation          : dU/dt + dF/dx + dG/dy = 0\n');
fprintf('Scheme            : finite volume, Rusanov flux, SSP-RK2, local time step\n');
fprintf('Reconstruction    : MUSCL minmod (order %d, first %d iterations first order)\n',order,nFirst);
fprintf('Grid              : %d x %d (%d cells)\n',Nx,Ny,Nx*Ny);
fprintf('Domain            : x [%.1f, %.1f]  y [%.1f, %.1f]\n',x_min,x_max,y_min,y_max);
fprintf('Mach number       : %.3f\n',M_inf);
fprintf('Angle of attack   : %.2f deg\n',aoa_deg);
fprintf('CFL               : %.2f\n',CFL);
fprintf('Maximum iterations: %d\n',maxIter);
fprintf('\n');

resHist = nan(maxIter,1);
ClHist  = nan(maxIter,1);
CdHist  = nan(maxIter,1);

res0 = [];
tic;

%% Time marching
for iter = 1:maxIter

    if iter <= nFirst
        ord = 1;
    else
        ord = order;
    end

    CFLn = CFL*min(1, 0.2 + 0.8*iter/rampIter);

    Uc0 = Uc;

    for stage = 1:2

        % Primitive variables
        rho = max(Uc(:,:,1), rho_floor);
        u   = Uc(:,:,2)./rho;
        v   = Uc(:,:,3)./rho;
        p   = max((gamma-1)*(Uc(:,:,4) - 0.5*rho.*(u.^2 + v.^2)), p_floor);

        Wp(2:Ny+1,2:Nx+1,1) = rho;
        Wp(2:Ny+1,2:Nx+1,2) = u;
        Wp(2:Ny+1,2:Nx+1,3) = v;
        Wp(2:Ny+1,2:Nx+1,4) = p;

        % Boundary conditions (ghost cells)
        Wp(:,1,:)    = repmat(reshape([rho_inf U_inf 0 p_inf],1,1,4),Ny+2,1,1);   % inlet
        Wp(1,:,:)    = repmat(reshape([rho_inf U_inf 0 p_inf],1,1,4),1,Nx+2,1);   % bottom
        Wp(Ny+2,:,:) = repmat(reshape([rho_inf U_inf 0 p_inf],1,1,4),1,Nx+2,1);   % top
        Wp(:,Nx+2,:) = Wp(:,Nx+1,:);                                               % outlet
        Wp(:,Nx+2,4) = p_inf;                                                      % outlet pressure

        % Local time step (frozen at the first stage)
        if stage == 1
            a = sqrt(gamma*p./rho);
            dtl = CFLn./((abs(u)+a)./DX + (abs(v)+a)./DY);
            dtl(solid) = 0;
        end

        Res = zeros(Ny,Nx,4);

        % Sweeps: dirn = 1 -> x faces, dirn = 2 -> y faces (arrays transposed)
        for dirn = 1:2

            if dirn == 1
                Wd = Wp;
                Sd = Sp;
                hp = dxp;
                cc = xcp;
            else
                Wd = permute(Wp,[2 1 3]);
                Wd = Wd(:,:,[1 3 2 4]);        % swap u <-> v: channel 2 is now normal velocity
                Sd = Sp.';
                hp = dyp;
                cc = ycp;
            end

            M2 = size(Wd,1);
            Wr = Wd(2:M2-1,:,:);               % transverse interior rows only
            Sr = Sd(2:M2-1,:);

            % MUSCL reconstruction with minmod limiter
            if ord == 2
                dq = diff(Wr,1,2)./diff(cc);
                sl = zeros(size(Wr));
                a1 = dq(:,1:end-1,:);
                b1 = dq(:,2:end,:);
                sl(:,2:end-1,:) = 0.5*(sign(a1)+sign(b1)).*min(abs(a1),abs(b1));

                % no reconstruction next to solid cells
                nbS = Sr(:,1:end-2) | Sr(:,2:end-1) | Sr(:,3:end);
                sl(:,2:end-1,:) = sl(:,2:end-1,:).*double(~nbS);

                WL = Wr(:,1:end-1,:) + 0.5*hp(1:end-1).*sl(:,1:end-1,:);
                WR = Wr(:,2:end,:)   - 0.5*hp(2:end).*sl(:,2:end,:);

                % fall back to first order where reconstruction is not admissible
                bad = WL(:,:,1) <= 0 | WL(:,:,4) <= 0 | WR(:,:,1) <= 0 | WR(:,:,4) <= 0;
                if any(bad(:))
                    for k = 1:4
                        c1 = Wr(:,1:end-1,k);
                        c2 = Wr(:,2:end,k);
                        tL = WL(:,:,k);  tR = WR(:,:,k);
                        tL(bad) = c1(bad);
                        tR(bad) = c2(bad);
                        WL(:,:,k) = tL;  WR(:,:,k) = tR;
                    end
                end
            else
                WL = Wr(:,1:end-1,:);
                WR = Wr(:,2:end,:);
            end

            % Solid faces: slip wall by mirroring the normal velocity of the fluid cell
            Wc1 = Wr(:,1:end-1,:);
            Wc2 = Wr(:,2:end,:);
            sL = Sr(:,1:end-1);
            sR = Sr(:,2:end);
            wallL = sL & ~sR;
            wallR = sR & ~sL;

            for k = 1:4
                tL = WL(:,:,k);  tR = WR(:,:,k);
                c1 = Wc1(:,:,k); c2 = Wc2(:,:,k);
                if k == 2
                    tL(wallL) = -c2(wallL);  tR(wallL) = c2(wallL);
                    tR(wallR) = -c1(wallR);  tL(wallR) = c1(wallR);
                else
                    tL(wallL) = c2(wallL);   tR(wallL) = c2(wallL);
                    tR(wallR) = c1(wallR);   tL(wallR) = c1(wallR);
                end
                WL(:,:,k) = tL;  WR(:,:,k) = tR;
            end

            % Rusanov flux
            rL = WL(:,:,1);  unL = WL(:,:,2);  utL = WL(:,:,3);  pL = WL(:,:,4);
            rR = WR(:,:,1);  unR = WR(:,:,2);  utR = WR(:,:,3);  pR = WR(:,:,4);

            EL = pL/(gamma-1) + 0.5*rL.*(unL.^2 + utL.^2);
            ER = pR/(gamma-1) + 0.5*rR.*(unR.^2 + utR.^2);

            UL = cat(3, rL, rL.*unL, rL.*utL, EL);
            UR = cat(3, rR, rR.*unR, rR.*utR, ER);

            FL = cat(3, rL.*unL, rL.*unL.^2 + pL, rL.*unL.*utL, unL.*(EL + pL));
            FR = cat(3, rR.*unR, rR.*unR.^2 + pR, rR.*unR.*utR, unR.*(ER + pR));

            smax = max(abs(unL) + sqrt(gamma*pL./rL), abs(unR) + sqrt(gamma*pR./rR));

            F = 0.5*(FL + FR) - 0.5*smax.*(UR - UL);

            % Flux divergence
            Rd = -(F(:,2:end,:) - F(:,1:end-1,:))./hp(2:end-1);

            if dirn == 1
                Res = Res + Rd;
            else
                Rd = Rd(:,:,[1 3 2 4]);        % swap momentum components back
                Res = Res + permute(Rd,[2 1 3]);
            end

        end

        % SSP-RK2 update
        if stage == 1
            Res1 = Res;
            Uc = Uc0 + dtl.*Res;
        else
            Uc = 0.5*Uc0 + 0.5*(Uc + dtl.*Res);
        end

    end

    % Residual (density, fluid cells only)
    r1 = Res1(:,:,1);
    res = sqrt(mean(r1(fluid).^2));
    if isempty(res0)
        res0 = res;
    end
    resHist(iter) = res/res0;

    % Divergence check
    if ~all(isfinite(Uc(:)))
        fprintf('\n');
        fprintf('SOLVER DIVERGED\n');
        fprintf('Iteration : %d\n',iter);
        fprintf('Time      : %.2f s\n',toc);
        break;
    end

    % Pressure force on the body (staircase surface, gauge pressure)
    rho = max(Uc(:,:,1), rho_floor);
    u   = Uc(:,:,2)./rho;
    v   = Uc(:,:,3)./rho;
    p   = max((gamma-1)*(Uc(:,:,4) - 0.5*rho.*(u.^2 + v.^2)), p_floor);

    Fx = sum(sum((p - p_inf).*wFx));
    Fy = sum(sum((p - p_inf).*wFy));
    CdHist(iter) = Fx/(q_inf*Lb);
    ClHist(iter) = Fy/(q_inf*Lb);

    % Live progress
    if mod(iter,printEvery) == 0
        fprintf(['Iteration: %5d / %5d   ' ...
                 'Progress: %6.2f %%   ' ...
                 'Residual: %.3e   ' ...
                 'Cl: %8.4f   Cd: %8.4f   ' ...
                 'Time: %.1f s\n'], ...
                 iter,maxIter,100*iter/maxIter,resHist(iter),ClHist(iter),CdHist(iter),toc);
    end

    % Convergence check
    if resHist(iter) < tol
        fprintf('\n');
        fprintf('Converged\n');
        fprintf('Iterations : %d\n',iter);
        fprintf('Residual   : %.3e\n',resHist(iter));
        fprintf('Time       : %.2f s\n',toc);
        break;
    end

end

resHist = resHist(1:iter);
ClHist  = ClHist(1:iter);
CdHist  = CdHist(1:iter);

if iter == maxIter && resHist(iter) >= tol
    fprintf('\n');
    fprintf('Maximum iterations reached\n');
    fprintf('Residual: %.3e\n',resHist(iter));
    fprintf('Time    : %.2f s\n',toc);
end

%% Post-processing
fprintf('\nCalculating flow fields...\n');

V    = sqrt(u.^2 + v.^2);
T    = p./(rho*R);
Mach = V./sqrt(gamma*R*T);
Cp   = (p - p_inf)/q_inf;

u(solid) = NaN;   v(solid) = NaN;   V(solid) = NaN;
p(solid) = NaN;   rho(solid) = NaN; T(solid) = NaN;
Mach(solid) = NaN;  Cp(solid) = NaN;

fprintf('\n');
fprintf('Solution complete\n');
fprintf('Angle of attack     : %.2f deg\n',aoa_deg);
fprintf('Freestream velocity : %.3f m/s\n',U_inf);
fprintf('Maximum velocity    : %.3f m/s\n',max(V(~isnan(V))));
fprintf('Minimum pressure    : %.2f Pa\n',min(p(~isnan(p))));
fprintf('Maximum pressure    : %.2f Pa\n',max(p(~isnan(p))));
fprintf('Maximum Mach number : %.4f\n',max(Mach(~isnan(Mach))));
fprintf('Lift coefficient    : %.4f  (based on Lb = %.2f)\n',ClHist(end),Lb);
fprintf('Drag coefficient    : %.4f\n',CdHist(end));

%% Plot window around the body
xlo = xc_b - 3;  xhi = xc_b + 6;
ylo = -3;        yhi = 3;

%% Plot 1: velocity magnitude
figure;
contourf(X,Y,V,30,'LineColor','none');
hold on;
contour(X,Y,double(solid),[0.5 0.5],'k','LineWidth',2);
axis equal;
xlim([xlo xhi]); ylim([ylo yhi]);
colorbar;
xlabel('x'); ylabel('y');
title(sprintf('Velocity Magnitude (AoA = %.1f deg)',aoa_deg));

%% Plot 2: pressure
figure;
contourf(X,Y,p,30,'LineColor','none');
hold on;
contour(X,Y,double(solid),[0.5 0.5],'k','LineWidth',2);
axis equal;
xlim([xlo xhi]); ylim([ylo yhi]);
colorbar;
xlabel('x'); ylabel('y');
title('Pressure');

%% Plot 3: Mach number
figure;
contourf(X,Y,Mach,30,'LineColor','none');
hold on;
contour(X,Y,double(solid),[0.5 0.5],'k','LineWidth',2);
axis equal;
xlim([xlo xhi]); ylim([ylo yhi]);
colorbar;
xlabel('x'); ylabel('y');
title('Mach Number');

%% Plot 4: velocity field
figure;
ix = find(x >= xlo & x <= xhi);
iy = find(y >= ylo & y <= yhi);
skip = 3;
ix = ix(1:skip:end);
iy = iy(1:skip:end);
quiver(X(iy,ix),Y(iy,ix),u(iy,ix),v(iy,ix));
hold on;
contour(X,Y,double(solid),[0.5 0.5],'k','LineWidth',2);
axis equal;
xlim([xlo xhi]); ylim([ylo yhi]);
xlabel('x'); ylabel('y');
title('Velocity Field');

%% Plot 5: convergence and force history
figure;
subplot(2,1,1);
semilogy(1:iter,resHist,'LineWidth',1.5);
grid on;
xlabel('Iteration'); ylabel('Relative density residual');
title('Convergence History');
subplot(2,1,2);
plot(1:iter,ClHist,'LineWidth',1.5); hold on;
plot(1:iter,CdHist,'LineWidth',1.5);
grid on;
xlabel('Iteration'); ylabel('Coefficient');
legend('C_l','C_d');

fprintf('\nAll plots generated.\n');
