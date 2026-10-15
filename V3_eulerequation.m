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

%% Body
shape = 'naca0012';          % 'rectangle', 'cylinder', 'ellipse', 'naca0012'

Rc   = 0.5;                  % cylinder radius
Lb   = 1.0;                  % rectangle / ellipse length, NACA chord
Hb   = 0.7;                  % rectangle / ellipse thickness
xc_b = 1.3;                  % body centre x (rotation centre)
yc_b = 0.0;                  % body centre y

switch lower(shape)
    case 'rectangle'
        Lref = Lb;           halfExt = 0.5*sqrt(Lb^2 + Hb^2);
    case 'cylinder'
        Lref = 2*Rc;         halfExt = Rc;
    case 'ellipse'
        Lref = Lb;           halfExt = 0.5*max(Lb,Hb);
    case 'naca0012'
        Lref = Lb;           halfExt = 0.5*Lb;
    otherwise
        error('Unknown shape: %s',shape);
end

%% Numerical methods
mode       = 'transient';    % 'steady' (local time step) or 'transient' (global time step)
flux       = 'hllc';         % 'rusanov', 'hll', 'hllc'
order      = 2;              % 1 = first order, 2 = MUSCL
limiter    = 'vanleer';      % 'minmod', 'vanleer', 'mc', 'vanalbada'
timeScheme = 'rk3';          % 'rk2', 'rk3'
farfield   = 'characteristic';   % 'characteristic' or 'freestream'

CFL        = 0.8;
nFirst     = 200;            % first-order start-up iterations
rampIter   = 300;            % CFL ramp length
maxIter    = 4000;
tol        = 1e-5;           % steady only: relative density residual
printEvery = 1;

rho_floor = 0.01*rho_inf;
p_floor   = 0.01*p_inf;

%% Transient / vortex development settings
tEnd_conv  = 50;             % end time in convective units  t* = t*U_inf/Lref
plotEvery  = 5;            % plot every N iterations
plotField  = 'vorticity';    % 'mach','pressure','density','velocity','vorticity','schlieren','entropy'
saveFrames = false;          % true = save frame_XXXXX.png
perturbAmp = 1e-2;           % start-up wake perturbation (fraction of U_inf), breaks symmetry
                             % set 0 for a perfectly symmetric start

transient = strcmpi(mode,'transient');
tEnd = tEnd_conv*Lref/U_inf;

%% Domain size (in units of Lref)
x_min = xc_b - 15*Lref;
x_max = xc_b + 30*Lref;      % long downstream length for the wake
y_min = yc_b - 15*Lref;
y_max = yc_b + 15*Lref;

%% Grid: uniform refined box (extended downstream for the wake) + geometric stretching
cellsPerLref = 20;           % resolution in the refined region
h0    = Lref/cellsPerLref;
r     = 1.05;                % stretching ratio outside the refined box
aUp   = 1.5*Lref;            % refined box upstream of body centre
aDown = 8*Lref;              % refined box downstream of body centre (vortex street)
aY    = 1.5*Lref;            % refined box half-height

if aY < halfExt || aUp < halfExt
    error('Refined box is smaller than the body');
end

% x direction
nx   = round((aUp+aDown)/h0);
xin  = linspace(xc_b-aUp, xc_b+aDown, nx+1);
LL   = (xc_b-aUp) - x_min;
LR   = x_max - (xc_b+aDown);
nL   = ceil(log(1 + LL*(r-1)/h0)/log(r));
nR   = ceil(log(1 + LR*(r-1)/h0)/log(r));
wL   = h0*r.^(1:nL);  wL = wL*(LL/sum(wL));
wR   = h0*r.^(1:nR);  wR = wR*(LR/sum(wR));
xf   = [fliplr((xc_b-aUp) - cumsum(wL)), xin, (xc_b+aDown) + cumsum(wR)];

% y direction
ny   = round(2*aY/h0);
yin  = linspace(yc_b-aY, yc_b+aY, ny+1);
LB   = (yc_b-aY) - y_min;
LT   = y_max - (yc_b+aY);
nB   = ceil(log(1 + LB*(r-1)/h0)/log(r));
nT   = ceil(log(1 + LT*(r-1)/h0)/log(r));
wB   = h0*r.^(1:nB);  wB = wB*(LB/sum(wB));
wT   = h0*r.^(1:nT);  wT = wT*(LT/sum(wT));
yf   = [fliplr((yc_b-aY) - cumsum(wB)), yin, (yc_b+aY) + cumsum(wT)];

x  = 0.5*(xf(1:end-1) + xf(2:end));
y  = 0.5*(yf(1:end-1) + yf(2:end));
dx = diff(xf);
dy = diff(yf);

Nx = numel(x);
Ny = numel(y);

[X,Y] = meshgrid(x,y);
DX = ones(Ny,1)*dx;
DY = dy(:)*ones(1,Nx);

dxp = [dx(1), dx, dx(end)];
dyp = [dy(1), dy, dy(end)];
xcp = [x(1)-dx(1), x, x(end)+dx(end)];
ycp = [y(1)-dy(1), y, y(end)+dy(end)];

%% Body mask (rotated clockwise by aoa, so positive aoa = nose up)
ca = cos(aoa);
sa = sin(aoa);

xb =  (X-xc_b)*ca - (Y-yc_b)*sa;
yb =  (X-xc_b)*sa + (Y-yc_b)*ca;

switch lower(shape)
    case 'rectangle'
        solid = (abs(xb) <= Lb/2 & abs(yb) <= Hb/2);
    case 'cylinder'
        solid = ((X-xc_b).^2 + (Y-yc_b).^2 <= Rc^2);
    case 'ellipse'
        solid = ((xb/(Lb/2)).^2 + (yb/(Hb/2)).^2 <= 1);
    case 'naca0012'
        xn = xb/Lb + 0.5;                       % 0 = leading edge, 1 = trailing edge
        tt = 0.12;
        yt = 5*tt*Lb*(0.2969*sqrt(max(xn,0)) - 0.1260*xn - 0.3516*xn.^2 ...
                      + 0.2843*xn.^3 - 0.1015*xn.^4);
        solid = (xn >= 0 & xn <= 1 & abs(yb) <= yt);
end
fluid = ~solid;

% Plot window (only these cells are drawn, which keeps the graphics light)
iw = find(x >= xc_b-3*Lref & x <= xc_b+12*Lref);
jw = find(y >= yc_b-3*Lref & y <= yc_b+3*Lref);

Sp = false(Ny+2,Nx+2);
Sp(2:Ny+1,2:Nx+1) = solid;

% Weights for the pressure force on the (staircase) body surface
sE = double(Sp(2:Ny+1,3:Nx+2));
sW = double(Sp(2:Ny+1,1:Nx));
sN = double(Sp(3:Ny+2,2:Nx+1));
sS = double(Sp(1:Ny,2:Nx+1));
wFx = double(fluid).*(sE - sW).*DY;
wFy = double(fluid).*(sN - sS).*DX;

%% Initial condition: uniform freestream (+ small wake perturbation)
Uc = zeros(Ny,Nx,4);
Uc(:,:,1) = rho_inf;
Uc(:,:,2) = rho_inf*U_inf;
Uc(:,:,3) = 0;
Uc(:,:,4) = p_inf/(gamma-1) + 0.5*rho_inf*U_inf^2;

if perturbAmp > 0
    vp = perturbAmp*U_inf*exp(-((X-(xc_b+1.5*Lref)).^2 + (Y-yc_b-0.25*Lref).^2)/(0.5*Lref)^2);
    vp(solid) = 0;
    Uc(:,:,3) = Uc(:,:,3) + rho_inf*vp;
end

Wp = zeros(Ny+2,Nx+2,4);     % padded primitive variables [rho u v p]
Wp(:,:,1) = rho_inf;
Wp(:,:,2) = U_inf;
Wp(:,:,3) = 0;
Wp(:,:,4) = p_inf;

% Wake probe
[~,ixp] = min(abs(x - (xc_b+1.5*Lref)));
[~,iyp] = min(abs(y - (yc_b+0.25*Lref)));

if strcmpi(timeScheme,'rk3')
    nStages = 3;
else
    nStages = 2;
end

%% Solver start
fprintf('\n');
fprintf('2D Euler Solver (Version 3)\n');
fprintf('Body              : %s (Lref = %.3f)\n',shape,Lref);
fprintf('Mode              : %s\n',mode);
fprintf('Riemann solver    : %s\n',flux);
fprintf('Reconstruction    : order %d, limiter %s (first %d iterations first order)\n',order,limiter,nFirst);
fprintf('Time scheme       : SSP-%s, CFL %.2f\n',upper(timeScheme),CFL);
fprintf('Far field         : %s\n',farfield);
fprintf('Grid              : %d x %d (%d cells), h0 = %.4f\n',Nx,Ny,Nx*Ny,h0);
fprintf('Domain            : x [%.1f, %.1f]  y [%.1f, %.1f]\n',x_min,x_max,y_min,y_max);
fprintf('Mach number       : %.3f\n',M_inf);
fprintf('Angle of attack   : %.2f deg\n',aoa_deg);
if transient
    fprintf('End time          : t* = %.1f, plot every %d iterations\n',tEnd_conv,plotEvery);
end
fprintf('\n');

resHist = nan(maxIter,1);
ClHist  = nan(maxIter,1);
CdHist  = nan(maxIter,1);
tHist   = nan(maxIter,1);
vpHist  = nan(maxIter,1);

res0 = [];
t = 0;
dtmin = 0;

if transient
    hf = figure('Name','Transient flow');
end

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

    for stage = 1:nStages

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
        if strcmpi(farfield,'characteristic')
            for side = 1:4
                switch side
                    case 1      % inlet (left)
                        nxb = -1;  nyb = 0;
                        bI_r = Wp(2:Ny+1,2,1);     bI_u = Wp(2:Ny+1,2,2);
                        bI_v = Wp(2:Ny+1,2,3);     bI_p = Wp(2:Ny+1,2,4);
                    case 2      % outlet (right)
                        nxb =  1;  nyb = 0;
                        bI_r = Wp(2:Ny+1,Nx+1,1);  bI_u = Wp(2:Ny+1,Nx+1,2);
                        bI_v = Wp(2:Ny+1,Nx+1,3);  bI_p = Wp(2:Ny+1,Nx+1,4);
                    case 3      % bottom
                        nxb = 0;   nyb = -1;
                        bI_r = Wp(2,2:Nx+1,1);     bI_u = Wp(2,2:Nx+1,2);
                        bI_v = Wp(2,2:Nx+1,3);     bI_p = Wp(2,2:Nx+1,4);
                    case 4      % top
                        nxb = 0;   nyb = 1;
                        bI_r = Wp(Ny+1,2:Nx+1,1);  bI_u = Wp(Ny+1,2:Nx+1,2);
                        bI_v = Wp(Ny+1,2:Nx+1,3);  bI_p = Wp(Ny+1,2:Nx+1,4);
                end

                aI  = sqrt(gamma*bI_p./bI_r);
                unI = bI_u*nxb + bI_v*nyb;
                unF = U_inf*nxb;

                Rp  = unI + 2*aI/(gamma-1);                 % outgoing invariant (interior)
                Rm  = unF - 2*a_inf/(gamma-1);              % incoming invariant (freestream)
                unB = 0.5*(Rp + Rm);
                aB  = max(0.25*(gamma-1)*(Rp - Rm), 0.1*a_inf);

                isOut = double(unB > 0);
                sI = bI_p./bI_r.^gamma;
                sF = p_inf/rho_inf^gamma;
                sB = isOut.*sI + (1-isOut)*sF;
                utx = isOut.*(bI_u - unI*nxb) + (1-isOut)*(U_inf - unF*nxb);
                uty = isOut.*(bI_v - unI*nyb) + (1-isOut)*(0     - unF*nyb);

                rB = (aB.^2./(gamma*sB)).^(1/(gamma-1));
                pB = rB.*aB.^2/gamma;
                uB = unB*nxb + utx;
                vB = unB*nyb + uty;

                switch side
                    case 1
                        Wp(2:Ny+1,1,1) = rB;  Wp(2:Ny+1,1,2) = uB;  Wp(2:Ny+1,1,3) = vB;  Wp(2:Ny+1,1,4) = pB;
                    case 2
                        Wp(2:Ny+1,Nx+2,1) = rB;  Wp(2:Ny+1,Nx+2,2) = uB;  Wp(2:Ny+1,Nx+2,3) = vB;  Wp(2:Ny+1,Nx+2,4) = pB;
                    case 3
                        Wp(1,2:Nx+1,1) = rB;  Wp(1,2:Nx+1,2) = uB;  Wp(1,2:Nx+1,3) = vB;  Wp(1,2:Nx+1,4) = pB;
                    case 4
                        Wp(Ny+2,2:Nx+1,1) = rB;  Wp(Ny+2,2:Nx+1,2) = uB;  Wp(Ny+2,2:Nx+1,3) = vB;  Wp(Ny+2,2:Nx+1,4) = pB;
                end
            end
        else
            Wp(:,1,:)    = repmat(reshape([rho_inf U_inf 0 p_inf],1,1,4),Ny+2,1,1);
            Wp(1,:,:)    = repmat(reshape([rho_inf U_inf 0 p_inf],1,1,4),1,Nx+2,1);
            Wp(Ny+2,:,:) = repmat(reshape([rho_inf U_inf 0 p_inf],1,1,4),1,Nx+2,1);
            Wp(:,Nx+2,:) = Wp(:,Nx+1,:);
            Wp(:,Nx+2,4) = p_inf;
        end

        % Time step (frozen at the first stage)
        if stage == 1
            a = sqrt(gamma*p./rho);
            dtl = CFLn./((abs(u)+a)./DX + (abs(v)+a)./DY);
            dtl(solid) = 0;
            if transient
                dtmin = min(dtl(fluid));            % global time step
                dtl = dtmin*double(fluid);
            end
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
                Wd = Wd(:,:,[1 3 2 4]);        % channel 2 becomes the normal velocity
                Sd = Sp.';
                hp = dyp;
                cc = ycp;
            end

            M2 = size(Wd,1);
            Wr = Wd(2:M2-1,:,:);
            Sr = Sd(2:M2-1,:);

            % MUSCL reconstruction
            if ord == 2
                dq = diff(Wr,1,2)./diff(cc);
                sl = zeros(size(Wr));
                a1 = dq(:,1:end-1,:);
                b1 = dq(:,2:end,:);

                switch lower(limiter)
                    case 'minmod'
                        phi = 0.5*(sign(a1)+sign(b1)).*min(abs(a1),abs(b1));
                    case 'vanleer'
                        phi = (a1.*abs(b1) + abs(a1).*b1)./(abs(a1) + abs(b1) + 1e-300);
                    case 'mc'
                        phi = 0.5*(sign(a1)+sign(b1)).*min(min(2*abs(a1),2*abs(b1)),0.5*abs(a1+b1));
                    case 'vanalbada'
                        phi = double(a1.*b1 > 0).*a1.*b1.*(a1+b1)./(a1.^2 + b1.^2 + 1e-300);
                    otherwise
                        error('Unknown limiter: %s',limiter);
                end
                sl(:,2:end-1,:) = phi;

                % no reconstruction next to solid cells
                nbS = Sr(:,1:end-2) | Sr(:,2:end-1) | Sr(:,3:end);
                sl(:,2:end-1,:) = sl(:,2:end-1,:).*double(~nbS);

                WL = Wr(:,1:end-1,:) + 0.5*hp(1:end-1).*sl(:,1:end-1,:);
                WR = Wr(:,2:end,:)   - 0.5*hp(2:end).*sl(:,2:end,:);

                % first-order fallback where reconstruction is not admissible
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

            % Solid faces: slip wall, mirror the normal velocity of the fluid cell
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

            % Left / right states
            rL = WL(:,:,1);  unL = WL(:,:,2);  utL = WL(:,:,3);  pL = WL(:,:,4);
            rR = WR(:,:,1);  unR = WR(:,:,2);  utR = WR(:,:,3);  pR = WR(:,:,4);

            EL = pL/(gamma-1) + 0.5*rL.*(unL.^2 + utL.^2);
            ER = pR/(gamma-1) + 0.5*rR.*(unR.^2 + utR.^2);

            UL = cat(3, rL, rL.*unL, rL.*utL, EL);
            UR = cat(3, rR, rR.*unR, rR.*utR, ER);

            FL = cat(3, rL.*unL, rL.*unL.^2 + pL, rL.*unL.*utL, unL.*(EL + pL));
            FR = cat(3, rR.*unR, rR.*unR.^2 + pR, rR.*unR.*utR, unR.*(ER + pR));

            aL = sqrt(gamma*pL./rL);
            aR = sqrt(gamma*pR./rR);

            % Riemann solver
            if strcmpi(flux,'rusanov')

                smax = max(abs(unL) + aL, abs(unR) + aR);
                F = 0.5*(FL + FR) - 0.5*smax.*(UR - UL);

            else

                % Roe averages -> Einfeldt wave speed estimates
                sq1 = sqrt(rL);  sq2 = sqrt(rR);
                w1  = sq1./(sq1 + sq2);  w2 = 1 - w1;
                HL  = (EL + pL)./rL;     HR = (ER + pR)./rR;
                unT = w1.*unL + w2.*unR;
                utT = w1.*utL + w2.*utR;
                HT  = w1.*HL  + w2.*HR;
                aT  = sqrt(max((gamma-1)*(HT - 0.5*(unT.^2 + utT.^2)), 1e-12));

                SL = min(unL - aL, unT - aT);
                SR = max(unR + aR, unT + aT);

                if strcmpi(flux,'hll')

                    Fh = (SR.*FL - SL.*FR + SL.*SR.*(UR - UL))./(SR - SL);
                    m1 = double(SL >= 0);
                    m3 = double(SR <= 0);
                    m2 = 1 - m1 - m3;
                    F  = m1.*FL + m2.*Fh + m3.*FR;

                else   % HLLC

                    Ss = (pR - pL + rL.*unL.*(SL - unL) - rR.*unR.*(SR - unR)) ...
                         ./(rL.*(SL - unL) - rR.*(SR - unR));

                    kL = rL.*(SL - unL)./(SL - Ss);
                    kR = rR.*(SR - unR)./(SR - Ss);

                    UsL = cat(3, kL, kL.*Ss, kL.*utL, ...
                              kL.*(EL./rL + (Ss - unL).*(Ss + pL./(rL.*(SL - unL)))));
                    UsR = cat(3, kR, kR.*Ss, kR.*utR, ...
                              kR.*(ER./rR + (Ss - unR).*(Ss + pR./(rR.*(SR - unR)))));

                    FsL = FL + SL.*(UsL - UL);
                    FsR = FR + SR.*(UsR - UR);

                    m1 = double(SL >= 0);
                    m2 = double(SL <  0 & Ss >= 0);
                    m3 = double(Ss <  0 & SR >  0);
                    m4 = double(SR <= 0);

                    F = m1.*FL + m2.*FsL + m3.*FsR + m4.*FR;

                end
            end

            % Flux divergence
            Rd = -(F(:,2:end,:) - F(:,1:end-1,:))./hp(2:end-1);

            if dirn == 1
                Res = Res + Rd;
            else
                Rd = Rd(:,:,[1 3 2 4]);
                Res = Res + permute(Rd,[2 1 3]);
            end

        end

        % SSP Runge-Kutta update
        if stage == 1
            Res1 = Res;
            Uc = Uc0 + dtl.*Res;
        elseif stage == 2
            if nStages == 3
                Uc = 0.75*Uc0 + 0.25*(Uc + dtl.*Res);
            else
                Uc = 0.5*Uc0 + 0.5*(Uc + dtl.*Res);
            end
        else
            Uc = Uc0/3 + (2/3)*(Uc + dtl.*Res);
        end

    end

    if transient
        t = t + dtmin;
    end
    tHist(iter) = t;

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

    % Updated primitives, force and probe
    rho = max(Uc(:,:,1), rho_floor);
    u   = Uc(:,:,2)./rho;
    v   = Uc(:,:,3)./rho;
    p   = max((gamma-1)*(Uc(:,:,4) - 0.5*rho.*(u.^2 + v.^2)), p_floor);

    Fx = sum(sum((p - p_inf).*wFx));
    Fy = sum(sum((p - p_inf).*wFy));
    CdHist(iter) = Fx/(q_inf*Lref);
    ClHist(iter) = Fy/(q_inf*Lref);
    vpHist(iter) = v(iyp,ixp);

    % Live progress
    if mod(iter,printEvery) == 0
        fprintf(['Iteration: %6d / %6d   ' ...
                 'Residual: %.3e   ' ...
                 'Cl: %8.4f   Cd: %8.4f   ' ...
                 't*: %7.2f   ' ...
                 'Time: %.1f s\n'], ...
                 iter,maxIter,resHist(iter),ClHist(iter),CdHist(iter),t*U_inf/Lref,toc);
    end

    % Transient: plot the evolving flow
    if transient && (mod(iter,plotEvery) == 0 || iter == 1)

        Vm = sqrt(u.^2 + v.^2);
        switch lower(plotField)
            case 'mach'
                fld = Vm./sqrt(gamma*p./rho);   ttl = 'Mach number';
            case 'pressure'
                fld = p;                        ttl = 'Pressure [Pa]';
            case 'density'
                fld = rho;                      ttl = 'Density [kg/m^3]';
            case 'velocity'
                fld = Vm;                       ttl = 'Velocity magnitude [m/s]';
            case 'vorticity'
                [dvdx,~] = gradient(v,x,y);
                [~,dudy] = gradient(u,x,y);
                fld = (dvdx - dudy)*Lref/U_inf; ttl = 'Vorticity \omega L_{ref}/U_\infty';
            case 'schlieren'
                [drdx,drdy] = gradient(rho,x,y);
                gm = sqrt(drdx.^2 + drdy.^2);
                fld = exp(-15*gm/max(gm(fluid)));  ttl = 'Numerical schlieren';
            case 'entropy'
                fld = p./p_inf.*(rho_inf./rho).^gamma - 1;  ttl = 'Entropy error (s/s_\infty - 1)';
            otherwise
                error('Unknown plotField: %s',plotField);
        end
        fld(solid) = NaN;

        figure(hf); clf;
        pcolor(X(jw,iw),Y(jw,iw),fld(jw,iw)); shading flat;
        if strcmpi(plotField,'schlieren')
            colormap(gray);
        elseif strcmpi(plotField,'vorticity')
            fw = fld(jw,iw);  vmax = max(abs(fw(~isnan(fw))));
            caxis([-0.3*vmax 0.3*vmax]);
        end
        hold on;
        contour(X(jw,iw),Y(jw,iw),double(solid(jw,iw)),[0.5 0.5],'k','LineWidth',2);
        axis equal;
        xlim([xc_b-3*Lref xc_b+12*Lref]); ylim([yc_b-3*Lref yc_b+3*Lref]);
        colorbar;
        xlabel('x'); ylabel('y');
        title(sprintf('%s   t* = %.2f   iter = %d',ttl,t*U_inf/Lref,iter));
        drawnow;

        if saveFrames
            print(hf,sprintf('frame_%05d.png',iter),'-dpng','-r120');
        end
    end

    % Transient end time
    if transient && t >= tEnd
        fprintf('\n');
        fprintf('End time reached\n');
        fprintf('Iterations : %d\n',iter);
        fprintf('t*         : %.2f\n',t*U_inf/Lref);
        fprintf('Time       : %.2f s\n',toc);
        break;
    end

    % Convergence check (steady only)
    if ~transient && resHist(iter) < tol
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
tHist   = tHist(1:iter);
vpHist  = vpHist(1:iter);

if iter == maxIter && ((~transient && resHist(iter) >= tol) || (transient && t < tEnd))
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
dS   = p./p_inf.*(rho_inf./rho).^gamma - 1;           % entropy error

[dvdx,~] = gradient(v,x,y);
[~,dudy] = gradient(u,x,y);
vort = (dvdx - dudy)*Lref/U_inf;

dS_L2  = sqrt(mean(dS(fluid).^2));
dS_max = max(dS(fluid));

u(solid) = NaN;   v(solid) = NaN;   V(solid) = NaN;
p(solid) = NaN;   rho(solid) = NaN; T(solid) = NaN;
Mach(solid) = NaN;  Cp(solid) = NaN;  dS(solid) = NaN;  vort(solid) = NaN;

fprintf('\n');
fprintf('Solution complete\n');
fprintf('Body                : %s, AoA = %.2f deg\n',shape,aoa_deg);
fprintf('Freestream velocity : %.3f m/s\n',U_inf);
fprintf('Maximum velocity    : %.3f m/s\n',max(V(~isnan(V))));
fprintf('Minimum pressure    : %.2f Pa\n',min(p(~isnan(p))));
fprintf('Maximum pressure    : %.2f Pa\n',max(p(~isnan(p))));
fprintf('Maximum Mach number : %.4f\n',max(Mach(~isnan(Mach))));
fprintf('Lift coefficient    : %.4f  (based on Lref = %.3f)\n',ClHist(end),Lref);
fprintf('Drag coefficient    : %.4f\n',CdHist(end));
fprintf('Entropy error L2    : %.3e\n',dS_L2);
fprintf('Entropy error max   : %.3e\n',dS_max);

% Vortex shedding frequency (transient)
if transient
    i2 = find(tHist >= 0.4*tHist(end));
    if numel(i2) > 64
        sigs  = {ClHist, vpHist};
        names = {'Cl signal', 'wake probe v'};
        tu = linspace(tHist(i2(1)), tHist(end), numel(i2));
        fs = 1/(tu(2) - tu(1));
        n  = numel(tu);
        win = 0.5*(1 - cos(2*pi*(0:n-1)/(n-1)));
        for ks = 1:2
            sg = interp1(tHist(i2), sigs{ks}(i2), tu);
            sg = (sg(:).' - mean(sg)).*win;
            Sf = abs(fft(sg));
            [pk,im] = max(Sf(2:floor(n/2)));
            f  = im*fs/n;
            fprintf('Shedding (%s): f = %.2f Hz, St = f Lref/U_inf = %.3f\n',names{ks},f,f*Lref/U_inf);
        end
        fprintf('Mean Cl (last 60%%) : %.4f\n',mean(ClHist(i2)));
        fprintf('Mean Cd (last 60%%) : %.4f\n',mean(CdHist(i2)));
        fprintf('(Spectrum resolution is limited by the record length; run longer for a sharper St.)\n');
    else
        fprintf('Run too short for a shedding frequency estimate.\n');
    end
end

%% Plot window around the body
xlo = xc_b - 3*Lref;  xhi = xc_b + 12*Lref;
ylo = yc_b - 3*Lref;  yhi = yc_b + 3*Lref;

%% Plot 1: Mach number
figure;
pcolor(X(jw,iw),Y(jw,iw),Mach(jw,iw)); shading flat;
hold on;
contour(X(jw,iw),Y(jw,iw),double(solid(jw,iw)),[0.5 0.5],'k','LineWidth',2);
axis equal; xlim([xlo xhi]); ylim([ylo yhi]);
colorbar; xlabel('x'); ylabel('y');
title('Mach Number');

%% Plot 2: pressure coefficient
figure;
pcolor(X(jw,iw),Y(jw,iw),Cp(jw,iw)); shading flat;
hold on;
contour(X(jw,iw),Y(jw,iw),double(solid(jw,iw)),[0.5 0.5],'k','LineWidth',2);
axis equal; xlim([xlo xhi]); ylim([ylo yhi]);
colorbar; xlabel('x'); ylabel('y');
title('Pressure Coefficient');

%% Plot 3: entropy error (numerical dissipation)
figure;
pcolor(X(jw,iw),Y(jw,iw),dS(jw,iw)); shading flat;
hold on;
contour(X(jw,iw),Y(jw,iw),double(solid(jw,iw)),[0.5 0.5],'k','LineWidth',2);
axis equal; xlim([xlo xhi]); ylim([ylo yhi]);
colorbar; xlabel('x'); ylabel('y');
title('Entropy Error  s/s_\infty - 1');

%% Plot 4: vorticity
figure;
vw = vort(jw,iw);  vmax = max(abs(vw(~isnan(vw))));
pcolor(X(jw,iw),Y(jw,iw),vort(jw,iw)); shading flat;
caxis([-0.3*vmax 0.3*vmax]);
hold on;
contour(X(jw,iw),Y(jw,iw),double(solid(jw,iw)),[0.5 0.5],'k','LineWidth',2);
axis equal; xlim([xlo xhi]); ylim([ylo yhi]);
colorbar; xlabel('x'); ylabel('y');
title('Vorticity  \omega L_{ref}/U_\infty');

%% Plot 5: velocity field (interpolated to a uniform grid)
figure;
xq = linspace(xlo,xhi,70);
yq = linspace(ylo,yhi,35);
[Xq,Yq] = meshgrid(xq,yq);
uq = interp2(x,y,u,Xq,Yq);
vq = interp2(x,y,v,Xq,Yq);
quiver(Xq,Yq,uq,vq);
hold on;
contour(X(jw,iw),Y(jw,iw),double(solid(jw,iw)),[0.5 0.5],'k','LineWidth',2);
axis equal; xlim([xlo xhi]); ylim([ylo yhi]);
xlabel('x'); ylabel('y');
title('Velocity Field');

%% Plot 6: histories
figure;
if transient
    xh = tHist*U_inf/Lref;   xlab = 't^* = t U_\infty / L_{ref}';
else
    xh = (1:iter)';          xlab = 'Iteration';
end
subplot(3,1,1);
semilogy(xh,resHist,'LineWidth',1.5);
grid on; xlabel(xlab); ylabel('Rel. density residual');
title('Convergence History');
subplot(3,1,2);
plot(xh,ClHist,'LineWidth',1.5); hold on;
plot(xh,CdHist,'LineWidth',1.5);
grid on; xlabel(xlab); ylabel('Coefficient');
legend('C_l','C_d');
subplot(3,1,3);
plot(xh,vpHist,'LineWidth',1.5);
grid on; xlabel(xlab); ylabel('v at wake probe [m/s]');

fprintf('\nAll plots generated.\n');
