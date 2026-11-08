function ns_solver(P)
%NS_SOLVER  2D laminar compressible Navier-Stokes, finite volume, body-fitted O-grid.
%
%  - Convective flux : Rusanov or HLLC, MUSCL (index space) with a limiter
%  - Viscous flux    : Green-Gauss cell gradients, corrected face gradients
%  - Time stepping   : SSP-RK2 / SSP-RK3, local time step (steady) or global (transient)
%  - Wall            : no-slip, adiabatic or isothermal
%  - Far field       : characteristic (Riemann invariant) boundary condition
%  - Mesh            : mesh_<geometry>.m  (O-grid, i = around the body, j = wall -> far field)
%
%  Non-dimensionalisation: rho_inf = 1, a_inf = 1, T_inf = 1, p_inf = 1/gamma,
%  U_inf = Mach, lengths in the units of the mesh files.
%  Run through ns_input.m.

P = set_defaults(P);

gamma = P.gamma;
M_inf = P.Mach;

rho_inf = 1;
a_inf   = 1;
T_inf   = 1;
p_inf   = 1/gamma;
U_inf   = M_inf;
q_inf   = 0.5*rho_inf*U_inf^2;

transient = strcmpi(P.mode,'transient');

%% Mesh
meshFun = str2func(['mesh_' lower(P.geometry)]);
[Xn,Yn,Lref] = meshFun(P);

al = P.aoa_deg*pi/180;                 % rotate the geometry (positive = nose up)
X =  Xn*cos(al) + Yn*sin(al);
Y = -Xn*sin(al) + Yn*cos(al);

Ni = size(X,1) - 1;
Nj = size(X,2) - 1;

XC = 0.25*(X(1:Ni,1:Nj) + X(2:Ni+1,1:Nj) + X(1:Ni,2:Nj+1) + X(2:Ni+1,2:Nj+1));
YC = 0.25*(Y(1:Ni,1:Nj) + Y(2:Ni+1,1:Nj) + Y(1:Ni,2:Nj+1) + Y(2:Ni+1,2:Nj+1));

Om = 0.5*abs((X(2:Ni+1,2:Nj+1) - X(1:Ni,1:Nj)).*(Y(1:Ni,2:Nj+1) - Y(2:Ni+1,1:Nj)) ...
           - (X(1:Ni,2:Nj+1) - X(2:Ni+1,1:Nj)).*(Y(2:Ni+1,2:Nj+1) - Y(1:Ni,1:Nj)));

if min(Om(:)) <= 1e-14
    error('Mesh has zero-area cells. Adjust Ni, Nj, h1_Lref or blendPower.');
end

% face area vectors (pointing toward increasing i and increasing j)
dxi = X(:,2:Nj+1) - X(:,1:Nj);
dyi = Y(:,2:Nj+1) - Y(:,1:Nj);
Six = -dyi;
Siy =  dxi;

dxj = X(2:Ni+1,:) - X(1:Ni,:);
dyj = Y(2:Ni+1,:) - Y(1:Ni,:);
Sjx =  dyj;
Sjy = -dxj;

C.gamma = gamma;
C.Pr    = P.Pr;
C.Ni    = Ni;
C.Nj    = Nj;
C.Om    = Om;
C.Six = Six;  C.Siy = Siy;
C.Sjx = Sjx;  C.Sjy = Sjy;
C.Ai  = hypot(Six,Siy);
C.Aj  = hypot(Sjx,Sjy);
C.nIx = Six./C.Ai;  C.nIy = Siy./C.Ai;
C.nJx = Sjx./C.Aj;  C.nJy = Sjy./C.Aj;

% padded cell centres (periodic in i, reflected ghost centres in j)
xcp = zeros(Ni+2,Nj+2);
ycp = zeros(Ni+2,Nj+2);
xcp(2:Ni+1,2:Nj+1) = XC;
ycp(2:Ni+1,2:Nj+1) = YC;
xw = 0.5*(X(1:Ni,1)    + X(2:Ni+1,1));     yw = 0.5*(Y(1:Ni,1)    + Y(2:Ni+1,1));
xf = 0.5*(X(1:Ni,Nj+1) + X(2:Ni+1,Nj+1));  yf = 0.5*(Y(1:Ni,Nj+1) + Y(2:Ni+1,Nj+1));
xcp(2:Ni+1,1)    = 2*xw - XC(:,1);   ycp(2:Ni+1,1)    = 2*yw - YC(:,1);
xcp(2:Ni+1,Nj+2) = 2*xf - XC(:,Nj);  ycp(2:Ni+1,Nj+2) = 2*yf - YC(:,Nj);
xcp(1,:) = xcp(Ni+1,:);  xcp(Ni+2,:) = xcp(2,:);
ycp(1,:) = ycp(Ni+1,:);  ycp(Ni+2,:) = ycp(2,:);

C.drIx = xcp(2:Ni+2,2:Nj+1) - xcp(1:Ni+1,2:Nj+1);
C.drIy = ycp(2:Ni+2,2:Nj+1) - ycp(1:Ni+1,2:Nj+1);
C.dsI  = hypot(C.drIx,C.drIy);
C.drJx = xcp(2:Ni+1,2:Nj+2) - xcp(2:Ni+1,1:Nj+1);
C.drJy = ycp(2:Ni+1,2:Nj+2) - ycp(2:Ni+1,1:Nj+1);
C.dsJ  = hypot(C.drJx,C.drJy);

% flow / fluid constants
C.mu0    = M_inf*Lref/P.Re;
C.suth   = strcmpi(P.viscosity,'sutherland');
C.S      = 110.4/P.Tinf_dim;
C.cpnd   = 1/(gamma-1);
C.limiter = P.limiter;
C.flux    = P.flux;
C.wall    = P.wall;
C.Tw      = P.Tw_ratio*T_inf;
C.rho_inf = rho_inf;
C.u_inf   = U_inf;
C.p_inf   = p_inf;
C.a_inf   = a_inf;
C.rho_floor = 0.01*rho_inf;
C.p_floor   = 0.01*p_inf;

nStages = 2;
if strcmpi(P.timeScheme,'rk3')
    nStages = 3;
end

tEnd = P.tEnd_conv*Lref/U_inf;

%% Initial condition
Uc = zeros(Ni,Nj,4);
Uc(:,:,1) = rho_inf;
Uc(:,:,2) = rho_inf*U_inf;
Uc(:,:,3) = 0;
Uc(:,:,4) = p_inf/(gamma-1) + 0.5*rho_inf*U_inf^2;

x0 = 1.5*Lref;                          % wake probe / perturbation position
y0 = 0.25*Lref;

if P.perturbAmp > 0
    vp = P.perturbAmp*U_inf*exp(-((XC-x0).^2 + (YC-y0).^2)/(0.5*Lref)^2);
    Uc(:,:,3) = Uc(:,:,3) + rho_inf*vp;
end

[~,ip] = min(reshape((XC-x0).^2 + (YC-y0).^2,[],1));
[ipi,ipj] = ind2sub([Ni Nj],ip);

%% Solver start
fprintf('\n');
fprintf('2D Navier-Stokes Solver (laminar, compressible)\n');
fprintf('Geometry          : %s (Lref = %.3f), AoA = %.2f deg\n',P.geometry,Lref,P.aoa_deg);
fprintf('Mesh              : %d x %d O-grid (%d cells), min cell area %.3e\n',Ni,Nj,Ni*Nj,min(Om(:)));
fprintf('Wall / far field  : %s wall, radius %.1f Lref\n',P.wall,P.Rfar_Lref);
fprintf('Mach / Reynolds   : %.3f / %.1f (Pr = %.2f, %s viscosity)\n',M_inf,P.Re,P.Pr,P.viscosity);
fprintf('Mode              : %s\n',P.mode);
fprintf('Convective flux   : %s, order %d, limiter %s\n',P.flux,P.order,P.limiter);
fprintf('Time scheme       : SSP-%s, CFL %.2f\n',upper(P.timeScheme),P.CFL);
if transient
    fprintf('End time          : t* = %.1f, plot every %d iterations\n',P.tEnd_conv,P.plotEvery);
end
fprintf('\n');

maxIter = P.maxIter;
resHist = nan(maxIter,1);
ClHist  = nan(maxIter,1);
CdHist  = nan(maxIter,1);
tHist   = nan(maxIter,1);
vpHist  = nan(maxIter,1);

res0  = [];
t     = 0;
dtmin = 0;
hf    = [];
if transient
    hf = figure('Name','Transient flow');
end

tic;

%% Time marching
for iter = 1:maxIter

    if iter <= P.nFirst
        ord = 1;
    else
        ord = P.order;
    end

    CFLn = P.CFL*min(1, 0.2 + 0.8*iter/P.rampIter);

    Uc0 = Uc;

    for stage = 1:nStages

        Wp = prim_from_cons(Uc,C);

        [Res,FvW] = residual(Wp,C,ord);

        if stage == 1
            dtl = local_dt(Wp,C,CFLn);
            if transient
                dtmin = min(dtl(:));
                dtl = dtmin*ones(Ni,Nj);
            end
            Res1 = Res;

            % forces on the body from the state at the start of the iteration
            pw  = Wp(2:Ni+1,2,4);
            Fxp = -sum(pw.*Sjx(:,1));
            Fyp = -sum(pw.*Sjy(:,1));
            Fxv =  sum(FvW(:,2));
            Fyv =  sum(FvW(:,3));
            CdHist(iter) = (Fxp + Fxv)/(q_inf*Lref);
            ClHist(iter) = (Fyp + Fyv)/(q_inf*Lref);
            CdP = Fxp/(q_inf*Lref);
            CdV = Fxv/(q_inf*Lref);
            vpHist(iter) = Wp(ipi+1,ipj+1,3);
        end

        if stage == 1
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

    r1  = Res1(:,:,1);
    res = sqrt(mean(r1(:).^2));
    if isempty(res0)
        res0 = res;
    end
    resHist(iter) = res/res0;

    if ~all(isfinite(Uc(:)))
        fprintf('\n');
        fprintf('SOLVER DIVERGED\n');
        fprintf('Iteration : %d\n',iter);
        fprintf('Time      : %.2f s\n',toc);
        break;
    end

    if mod(iter,P.printEvery) == 0
        fprintf(['Iteration: %6d / %6d   ' ...
                 'Residual: %.3e   ' ...
                 'Cl: %8.4f   Cd: %8.4f (p %.4f, v %.4f)   ' ...
                 't*: %7.2f   ' ...
                 'Time: %.1f s\n'], ...
                 iter,maxIter,resHist(iter),ClHist(iter),CdHist(iter),CdP,CdV,t*U_inf/Lref,toc);
    end

    % Transient: plot the evolving flow
    if transient && (mod(iter,P.plotEvery) == 0 || iter == 1)
        Wq = prim_from_cons(Uc,C);
        [fld,ttl] = plot_field(Wq,C,P.plotField,Lref,U_inf);
        figure(hf); clf;
        pcolor(X,Y,cell2node(fld)); shading flat;
        hold on;
        plot(X(:,1),Y(:,1),'k','LineWidth',1.5);
        axis equal;
        xlim([-3*Lref 12*Lref]); ylim([-4*Lref 4*Lref]);
        colorbar;
        if strcmpi(P.plotField,'vorticity')
            vm = 0.3*max(abs(fld(:)));
            caxis([-vm vm]);
        end
        xlabel('x'); ylabel('y');
        title(sprintf('%s   t* = %.2f   iter = %d',ttl,t*U_inf/Lref,iter));
        drawnow;
        if P.saveFrames
            print(hf,sprintf('frame_%05d.png',iter),'-dpng','-r120');
        end
    end

    if transient && t >= tEnd
        fprintf('\n');
        fprintf('End time reached\n');
        fprintf('Iterations : %d\n',iter);
        fprintf('t*         : %.2f\n',t*U_inf/Lref);
        fprintf('Time       : %.2f s\n',toc);
        break;
    end

    if ~transient && resHist(iter) < P.tol
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

if iter == maxIter && ((~transient && resHist(iter) >= P.tol) || (transient && t < tEnd))
    fprintf('\n');
    fprintf('Maximum iterations reached\n');
    fprintf('Residual: %.3e\n',resHist(iter));
    fprintf('Time    : %.2f s\n',toc);
end

%% Post-processing
fprintf('\nCalculating flow fields...\n');

Wp = prim_from_cons(Uc,C);
[Res,FvW] = residual(Wp,C,P.order);     %#ok<ASGLU>

rho = Wp(2:Ni+1,2:Nj+1,1);
u   = Wp(2:Ni+1,2:Nj+1,2);
v   = Wp(2:Ni+1,2:Nj+1,3);
p   = Wp(2:Ni+1,2:Nj+1,4);
Tn  = gamma*p./rho;
V   = hypot(u,v);
Mach = V./sqrt(Tn);                     % a = sqrt(T) in these units
Cp   = (p - p_inf)/q_inf;

[gvx,~]  = green_gauss(Wp(:,:,3),C);
[~,guy]  = green_gauss(Wp(:,:,2),C);
vort = (gvx(2:Ni+1,2:Nj+1) - guy(2:Ni+1,2:Nj+1))*Lref/U_inf;

% wall distributions
pw   = Wp(2:Ni+1,2,4);
Cp_w = (pw - p_inf)/q_inf;
tvec = [FvW(:,2)./C.Aj(:,1), FvW(:,3)./C.Aj(:,1)];
nxw  = C.nJx(:,1);  nyw = C.nJy(:,1);
tn   = tvec(:,1).*nxw + tvec(:,2).*nyw;
ttx  = tvec(:,1) - tn.*nxw;
tty  = tvec(:,2) - tn.*nyw;
Cf_w = (ttx.*(-nyw) + tty.*nxw)/q_inf;  % along the counter-clockwise direction
ds_w = C.Aj(:,1);
s_w  = (cumsum(ds_w) - 0.5*ds_w)/sum(ds_w);

fprintf('\n');
fprintf('Solution complete\n');
fprintf('Lift coefficient    : %.4f  (based on Lref = %.3f)\n',ClHist(end),Lref);
fprintf('Drag coefficient    : %.4f  (pressure %.4f, friction %.4f)\n',CdHist(end),CdP,CdV);
fprintf('Maximum Mach number : %.4f\n',max(Mach(:)));
fprintf('Minimum Cp          : %.4f\n',min(Cp_w));

% recirculation length (rear symmetry line, only for zero angle of attack)
if abs(P.aoa_deg) < 1e-9
    xl = 0.5*(XC(1,:) + XC(Ni,:));
    ul = 0.5*(u(1,:)  + u(Ni,:));
    k  = find(ul(1:end-1) < 0 & ul(2:end) >= 0,1);
    if ~isempty(k)
        xz = xl(k) - ul(k)*(xl(k+1) - xl(k))/(ul(k+1) - ul(k));
        fprintf('Recirculation length: %.3f Lref (from the rear of the body)\n',(xz - X(1,1))/Lref);
    else
        fprintf('Recirculation length: no reversed flow found on the rear symmetry line\n');
    end
end

% vortex shedding frequency (transient)
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
            [~,im] = max(Sf(2:floor(n/2)));
            f  = im*fs/n;
            fprintf('Shedding (%s): St = f Lref/U_inf = %.3f\n',names{ks},f*Lref/U_inf);
        end
        fprintf('Mean Cd (last 60%%) : %.4f\n',mean(CdHist(i2)));
        fprintf('(Spectrum resolution is limited by the record length; run longer for a sharper St.)\n');
    else
        fprintf('Run too short for a shedding frequency estimate.\n');
    end
end

%% Plots
xlo = -3*Lref;  xhi = 12*Lref;
ylo = -4*Lref;  yhi =  4*Lref;

figure;
pcolor(X,Y,cell2node(Mach)); shading flat; hold on;
plot(X(:,1),Y(:,1),'k','LineWidth',1.5);
axis equal; xlim([xlo xhi]); ylim([ylo yhi]); colorbar;
xlabel('x'); ylabel('y'); title('Mach Number');

figure;
pcolor(X,Y,cell2node(Cp)); shading flat; hold on;
plot(X(:,1),Y(:,1),'k','LineWidth',1.5);
axis equal; xlim([xlo xhi]); ylim([ylo yhi]); colorbar;
xlabel('x'); ylabel('y'); title('Pressure Coefficient');

figure;
pcolor(X,Y,cell2node(vort)); shading flat; hold on;
plot(X(:,1),Y(:,1),'k','LineWidth',1.5);
axis equal; xlim([xlo xhi]); ylim([ylo yhi]);
vm = 0.3*max(abs(vort(:)));
caxis([-vm vm]); colorbar;
xlabel('x'); ylabel('y'); title('Vorticity  \omega L_{ref}/U_\infty');

figure;
plot(X(:,1:3:end),Y(:,1:3:end),'k-','LineWidth',0.5); hold on;
plot(X(1:3:end,:).',Y(1:3:end,:).','k-','LineWidth',0.5);
axis equal; xlim([-1.5*Lref 2.5*Lref]); ylim([-1.5*Lref 1.5*Lref]);
xlabel('x'); ylabel('y'); title('Mesh (near the body)');

figure;
subplot(2,1,1);
plot(s_w,Cp_w,'LineWidth',1.5); grid on;
xlabel('surface position (rear point, counter-clockwise)'); ylabel('C_p');
title('Wall pressure and skin friction');
subplot(2,1,2);
plot(s_w,Cf_w,'LineWidth',1.5); grid on;
xlabel('surface position (rear point, counter-clockwise)'); ylabel('C_f');

figure;
if transient
    xh = tHist*U_inf/Lref;   xlab = 't^* = t U_\infty / L_{ref}';
else
    xh = (1:iter)';          xlab = 'Iteration';
end
subplot(3,1,1);
semilogy(xh,resHist,'LineWidth',1.5); grid on;
xlabel(xlab); ylabel('Rel. density residual'); title('Convergence History');
subplot(3,1,2);
plot(xh,ClHist,'LineWidth',1.5); hold on;
plot(xh,CdHist,'LineWidth',1.5); grid on;
xlabel(xlab); ylabel('Coefficient'); legend('C_l','C_d');
subplot(3,1,3);
plot(xh,vpHist,'LineWidth',1.5); grid on;
xlabel(xlab); ylabel('v at wake probe');

fprintf('\nAll plots generated.\n');

end


%% ========================================================================
%  Subfunctions
%  ========================================================================

function P = set_defaults(P)
D.geometry   = 'cylinder';
D.Rc         = 0.5;
D.Lb         = 1.0;
D.Hb         = 0.7;
D.naca_t     = 0.12;
D.aoa_deg    = 0;
D.Mach       = 0.3;
D.Re         = 40;
D.Pr         = 0.72;
D.gamma      = 1.4;
D.viscosity  = 'constant';
D.Tinf_dim   = 288.15;
D.wall       = 'adiabatic';
D.Tw_ratio   = 1.0;
D.Ni         = 160;
D.Nj         = 90;
D.h1_Lref    = 0.015;
D.Rfar_Lref  = 20;
D.blendPower = 1.5;
D.flux       = 'hllc';
D.order      = 2;
D.limiter    = 'vanleer';
D.timeScheme = 'rk3';
D.CFL        = 0.8;
D.nFirst     = 200;
D.rampIter   = 300;
D.mode       = 'steady';
D.maxIter    = 30000;
D.tol        = 1e-6;
D.tEnd_conv  = 100;
D.printEvery = 200;
D.plotEvery  = 500;
D.plotField  = 'vorticity';
D.saveFrames = false;
D.perturbAmp = 1e-2;
fn = fieldnames(D);
for k = 1:numel(fn)
    if ~isfield(P,fn{k})
        P.(fn{k}) = D.(fn{k});
    end
end
end


function Wp = prim_from_cons(Uc,C)
Ni = C.Ni;  Nj = C.Nj;  g = C.gamma;
rho = max(Uc(:,:,1),C.rho_floor);
u   = Uc(:,:,2)./rho;
v   = Uc(:,:,3)./rho;
p   = max((g-1)*(Uc(:,:,4) - 0.5*rho.*(u.^2 + v.^2)),C.p_floor);
Wp = zeros(Ni+2,Nj+2,4);
Wp(2:Ni+1,2:Nj+1,1) = rho;
Wp(2:Ni+1,2:Nj+1,2) = u;
Wp(2:Ni+1,2:Nj+1,3) = v;
Wp(2:Ni+1,2:Nj+1,4) = p;
Wp = apply_bc(Wp,C);
end


function Wp = apply_bc(Wp,C)
Ni = C.Ni;  Nj = C.Nj;  g = C.gamma;

% wall ghost cells (j = 1): no-slip, zero pressure gradient, adiabatic or isothermal
r1 = Wp(2:Ni+1,2,1);  u1 = Wp(2:Ni+1,2,2);
v1 = Wp(2:Ni+1,2,3);  p1 = Wp(2:Ni+1,2,4);
T1 = g*p1./r1;
if strcmpi(C.wall,'isothermal')
    Tg = max(2*C.Tw - T1, 0.2*C.Tw);
else
    Tg = T1;
end
Wp(2:Ni+1,1,1) = g*p1./Tg;
Wp(2:Ni+1,1,2) = -u1;
Wp(2:Ni+1,1,3) = -v1;
Wp(2:Ni+1,1,4) = p1;

% far-field ghost cells (j = Nj+2): characteristic boundary condition
nx = C.nJx(:,Nj+1);
ny = C.nJy(:,Nj+1);
rI = Wp(2:Ni+1,Nj+1,1);  uI = Wp(2:Ni+1,Nj+1,2);
vI = Wp(2:Ni+1,Nj+1,3);  pI = Wp(2:Ni+1,Nj+1,4);

aI  = sqrt(g*pI./rI);
unI = uI.*nx + vI.*ny;
unF = C.u_inf*nx;
Rp  = unI + 2*aI/(g-1);
Rm  = unF - 2*C.a_inf/(g-1);
unB = 0.5*(Rp + Rm);
aB  = max(0.25*(g-1)*(Rp - Rm), 0.1*C.a_inf);

isOut = double(unB > 0);
sI  = pI./rI.^g;
sF  = C.p_inf/C.rho_inf^g;
sB  = isOut.*sI + (1-isOut)*sF;
utx = isOut.*(uI - unI.*nx) + (1-isOut).*(C.u_inf - unF.*nx);
uty = isOut.*(vI - unI.*ny) + (1-isOut).*(0 - unF.*ny);

rB = (aB.^2./(g*sB)).^(1/(g-1));
pB = rB.*aB.^2/g;
Wp(2:Ni+1,Nj+2,1) = rB;
Wp(2:Ni+1,Nj+2,2) = unB.*nx + utx;
Wp(2:Ni+1,Nj+2,3) = unB.*ny + uty;
Wp(2:Ni+1,Nj+2,4) = pB;

% periodic in i
Wp(1,:,:)    = Wp(Ni+1,:,:);
Wp(Ni+2,:,:) = Wp(2,:,:);
end


function mu = viscosity(T,C)
if C.suth
    mu = C.mu0*T.^1.5.*(1 + C.S)./(T + C.S);
else
    mu = C.mu0*ones(size(T));
end
end


function dtl = local_dt(Wp,C,CFLn)
Ni = C.Ni;  Nj = C.Nj;  g = C.gamma;
rho = Wp(2:Ni+1,2:Nj+1,1);
u   = Wp(2:Ni+1,2:Nj+1,2);
v   = Wp(2:Ni+1,2:Nj+1,3);
p   = Wp(2:Ni+1,2:Nj+1,4);
a   = sqrt(g*p./rho);
mu  = viscosity(g*p./rho,C);

Sxi = 0.5*(C.Six(1:Ni,:) + C.Six(2:Ni+1,:));
Syi = 0.5*(C.Siy(1:Ni,:) + C.Siy(2:Ni+1,:));
Sxj = 0.5*(C.Sjx(:,1:Nj) + C.Sjx(:,2:Nj+1));
Syj = 0.5*(C.Sjy(:,1:Nj) + C.Sjy(:,2:Nj+1));

Lci = abs(u.*Sxi + v.*Syi) + a.*hypot(Sxi,Syi);
Lcj = abs(u.*Sxj + v.*Syj) + a.*hypot(Sxj,Syj);

nu  = max(4/3,g/C.Pr)*mu./rho;
Lvi = nu.*(Sxi.^2 + Syi.^2)./C.Om;
Lvj = nu.*(Sxj.^2 + Syj.^2)./C.Om;

dtl = CFLn*C.Om./(Lci + Lcj + 4*(Lvi + Lvj));
end


function [gx,gy] = green_gauss(phi,C)
% Green-Gauss cell gradients, returned on the padded array
Ni = C.Ni;  Nj = C.Nj;
phiI = 0.5*(phi(1:Ni+1,2:Nj+1) + phi(2:Ni+2,2:Nj+1));
phiJ = 0.5*(phi(2:Ni+1,1:Nj+1) + phi(2:Ni+1,2:Nj+2));

ax = (phiI(2:Ni+1,:).*C.Six(2:Ni+1,:) - phiI(1:Ni,:).*C.Six(1:Ni,:) ...
    + phiJ(:,2:Nj+1).*C.Sjx(:,2:Nj+1) - phiJ(:,1:Nj).*C.Sjx(:,1:Nj))./C.Om;
ay = (phiI(2:Ni+1,:).*C.Siy(2:Ni+1,:) - phiI(1:Ni,:).*C.Siy(1:Ni,:) ...
    + phiJ(:,2:Nj+1).*C.Sjy(:,2:Nj+1) - phiJ(:,1:Nj).*C.Sjy(:,1:Nj))./C.Om;

gx = zeros(Ni+2,Nj+2);
gy = zeros(Ni+2,Nj+2);
gx(2:Ni+1,2:Nj+1) = ax;
gy(2:Ni+1,2:Nj+1) = ay;
gx(2:Ni+1,1) = ax(:,1);   gx(2:Ni+1,Nj+2) = ax(:,Nj);
gy(2:Ni+1,1) = ay(:,1);   gy(2:Ni+1,Nj+2) = ay(:,Nj);
gx(1,:) = gx(Ni+1,:);     gx(Ni+2,:) = gx(2,:);
gy(1,:) = gy(Ni+1,:);     gy(Ni+2,:) = gy(2,:);
end


function [gx,gy] = face_grad(phiL,phiR,gLx,gLy,gRx,gRy,drx,dry,ds)
% face gradient: averaged cell gradients, corrected along the cell-to-cell direction
ax = 0.5*(gLx + gRx);
ay = 0.5*(gLy + gRy);
ex = drx./ds;
ey = dry./ds;
corr = (phiR - phiL)./ds - (ax.*ex + ay.*ey);
gx = ax + corr.*ex;
gy = ay + corr.*ey;
end


function Hv = visc_flux(uf,vf,muf,kf,ux,uy,vx,vy,Tx,Ty,Sx,Sy)
% viscous flux integrated over the face (area vector Sx, Sy)
div = ux + vy;
txx = muf.*(2*ux - (2/3)*div);
tyy = muf.*(2*vy - (2/3)*div);
txy = muf.*(uy + vx);
Hx = txx.*Sx + txy.*Sy;
Hy = txy.*Sx + tyy.*Sy;
He = (uf.*txx + vf.*txy + kf.*Tx).*Sx + (uf.*txy + vf.*tyy + kf.*Ty).*Sy;
Hv = cat(3,zeros(size(Hx)),Hx,Hy,He);
end


function [WL,WR] = recon(Wr,ord,limiter,periodic)
% MUSCL reconstruction along dimension 2 (index space)
if ord == 2
    dq = diff(Wr,1,2);
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
    sl = zeros(size(Wr));
    sl(:,2:end-1,:) = phi;
    if periodic
        sl(:,1,:)   = sl(:,end-1,:);
        sl(:,end,:) = sl(:,2,:);
    end
    WL = Wr(:,1:end-1,:) + 0.5*sl(:,1:end-1,:);
    WR = Wr(:,2:end,:)   - 0.5*sl(:,2:end,:);

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
end


function F = conv_flux(WL,WR,nx,ny,C)
% Rusanov or HLLC flux through faces with unit normal (nx,ny); returns flux per unit area
g = C.gamma;

rL = WL(:,:,1);  uL = WL(:,:,2);  vL = WL(:,:,3);  pL = WL(:,:,4);
rR = WR(:,:,1);  uR = WR(:,:,2);  vR = WR(:,:,3);  pR = WR(:,:,4);

unL =  uL.*nx + vL.*ny;   utL = -uL.*ny + vL.*nx;
unR =  uR.*nx + vR.*ny;   utR = -uR.*ny + vR.*nx;

EL = pL/(g-1) + 0.5*rL.*(unL.^2 + utL.^2);
ER = pR/(g-1) + 0.5*rR.*(unR.^2 + utR.^2);

UL = cat(3,rL,rL.*unL,rL.*utL,EL);
UR = cat(3,rR,rR.*unR,rR.*utR,ER);
FL = cat(3,rL.*unL,rL.*unL.^2 + pL,rL.*unL.*utL,unL.*(EL + pL));
FR = cat(3,rR.*unR,rR.*unR.^2 + pR,rR.*unR.*utR,unR.*(ER + pR));

aL = sqrt(g*pL./rL);
aR = sqrt(g*pR./rR);

if strcmpi(C.flux,'rusanov')

    smax = max(abs(unL) + aL,abs(unR) + aR);
    F = 0.5*(FL + FR) - 0.5*smax.*(UR - UL);

else

    sq1 = sqrt(rL);  sq2 = sqrt(rR);
    w1  = sq1./(sq1 + sq2);  w2 = 1 - w1;
    HL  = (EL + pL)./rL;     HR = (ER + pR)./rR;
    unT = w1.*unL + w2.*unR;
    utT = w1.*utL + w2.*utR;
    HT  = w1.*HL  + w2.*HR;
    aT  = sqrt(max((g-1)*(HT - 0.5*(unT.^2 + utT.^2)),1e-12));

    SL = min(unL - aL,unT - aT);
    SR = max(unR + aR,unT + aT);

    Ss = (pR - pL + rL.*unL.*(SL - unL) - rR.*unR.*(SR - unR)) ...
         ./(rL.*(SL - unL) - rR.*(SR - unR));

    kL = rL.*(SL - unL)./(SL - Ss);
    kR = rR.*(SR - unR)./(SR - Ss);

    UsL = cat(3,kL,kL.*Ss,kL.*utL, ...
              kL.*(EL./rL + (Ss - unL).*(Ss + pL./(rL.*(SL - unL)))));
    UsR = cat(3,kR,kR.*Ss,kR.*utR, ...
              kR.*(ER./rR + (Ss - unR).*(Ss + pR./(rR.*(SR - unR)))));

    FsL = FL + SL.*(UsL - UL);
    FsR = FR + SR.*(UsR - UR);

    m1 = double(SL >= 0);
    m2 = double(SL <  0 & Ss >= 0);
    m3 = double(Ss <  0 & SR >  0);
    m4 = double(SR <= 0);

    F = m1.*FL + m2.*FsL + m3.*FsR + m4.*FR;

end

% rotate the momentum flux back to the x-y frame
Fn = F(:,:,2);
Ft = F(:,:,3);
F  = cat(3,F(:,:,1),Fn.*nx - Ft.*ny,Fn.*ny + Ft.*nx,F(:,:,4));
end


function [R,FvW] = residual(Wp,C,ord)
Ni = C.Ni;  Nj = C.Nj;  g = C.gamma;

rho = Wp(:,:,1);  u = Wp(:,:,2);  v = Wp(:,:,3);  p = Wp(:,:,4);
T   = g*p./rho;
mu  = viscosity(T,C);
kk  = mu*C.cpnd/C.Pr;

[gux,guy] = green_gauss(u,C);
[gvx,gvy] = green_gauss(v,C);
[gTx,gTy] = green_gauss(T,C);

% ---- i-faces (between cells i-1 and i), size (Ni+1) x Nj
iL = 1:Ni+1;  iR = 2:Ni+2;  jj = 2:Nj+1;
[uxI,uyI] = face_grad(u(iL,jj),u(iR,jj),gux(iL,jj),guy(iL,jj),gux(iR,jj),guy(iR,jj),C.drIx,C.drIy,C.dsI);
[vxI,vyI] = face_grad(v(iL,jj),v(iR,jj),gvx(iL,jj),gvy(iL,jj),gvx(iR,jj),gvy(iR,jj),C.drIx,C.drIy,C.dsI);
[TxI,TyI] = face_grad(T(iL,jj),T(iR,jj),gTx(iL,jj),gTy(iL,jj),gTx(iR,jj),gTy(iR,jj),C.drIx,C.drIy,C.dsI);
HvI = visc_flux(0.5*(u(iL,jj)+u(iR,jj)),0.5*(v(iL,jj)+v(iR,jj)), ...
                0.5*(mu(iL,jj)+mu(iR,jj)),0.5*(kk(iL,jj)+kk(iR,jj)), ...
                uxI,uyI,vxI,vyI,TxI,TyI,C.Six,C.Siy);

% ---- j-faces (between cells j-1 and j), size Ni x (Nj+1)
ii = 2:Ni+1;  jL = 1:Nj+1;  jR = 2:Nj+2;
[uxJ,uyJ] = face_grad(u(ii,jL),u(ii,jR),gux(ii,jL),guy(ii,jL),gux(ii,jR),guy(ii,jR),C.drJx,C.drJy,C.dsJ);
[vxJ,vyJ] = face_grad(v(ii,jL),v(ii,jR),gvx(ii,jL),gvy(ii,jL),gvx(ii,jR),gvy(ii,jR),C.drJx,C.drJy,C.dsJ);
[TxJ,TyJ] = face_grad(T(ii,jL),T(ii,jR),gTx(ii,jL),gTy(ii,jL),gTx(ii,jR),gTy(ii,jR),C.drJx,C.drJy,C.dsJ);
HvJ = visc_flux(0.5*(u(ii,jL)+u(ii,jR)),0.5*(v(ii,jL)+v(ii,jR)), ...
                0.5*(mu(ii,jL)+mu(ii,jR)),0.5*(kk(ii,jL)+kk(ii,jR)), ...
                uxJ,uyJ,vxJ,vyJ,TxJ,TyJ,C.Sjx,C.Sjy);

FvW = squeeze(HvJ(:,1,:));              % viscous wall flux (Ni x 4)

% ---- convective flux, i-direction (arrays transposed so the sweep is along dim 2)
Wt = permute(Wp,[2 1 3]);
[WL,WR] = recon(Wt(2:Nj+1,:,:),ord,C.limiter,true);
Fc = conv_flux(WL,WR,C.nIx.',C.nIy.',C);
FcI = permute(Fc,[2 1 3]).*C.Ai;

% ---- convective flux, j-direction
[WL,WR] = recon(Wp(2:Ni+1,:,:),ord,C.limiter,false);
FcJ = conv_flux(WL,WR,C.nJx,C.nJy,C).*C.Aj;

% wall face: pressure flux only
pw = Wp(2:Ni+1,2,4);
FcJ(:,1,1) = 0;
FcJ(:,1,2) = pw.*C.Sjx(:,1);
FcJ(:,1,3) = pw.*C.Sjy(:,1);
FcJ(:,1,4) = 0;

HI = FcI - HvI;
HJ = FcJ - HvJ;

R = -(HI(2:Ni+1,:,:) - HI(1:Ni,:,:) + HJ(:,2:Nj+1,:) - HJ(:,1:Nj,:))./C.Om;
end


function Cn = cell2node(f)
% average cell values to the nodes (periodic in i, replicated in j)
[Ni,Nj] = size(f);
fp = zeros(Ni+2,Nj+2);
fp(2:Ni+1,2:Nj+1) = f;
fp(2:Ni+1,1)    = f(:,1);
fp(2:Ni+1,Nj+2) = f(:,Nj);
fp(1,:)    = fp(Ni+1,:);
fp(Ni+2,:) = fp(2,:);
Cn = 0.25*(fp(1:Ni+1,1:Nj+1) + fp(2:Ni+2,1:Nj+1) + fp(1:Ni+1,2:Nj+2) + fp(2:Ni+2,2:Nj+2));
end


function [fld,ttl] = plot_field(Wp,C,name,Lref,U_inf)
Ni = C.Ni;  Nj = C.Nj;  g = C.gamma;
rho = Wp(2:Ni+1,2:Nj+1,1);
u   = Wp(2:Ni+1,2:Nj+1,2);
v   = Wp(2:Ni+1,2:Nj+1,3);
p   = Wp(2:Ni+1,2:Nj+1,4);
switch lower(name)
    case 'mach'
        fld = hypot(u,v)./sqrt(g*p./rho);   ttl = 'Mach number';
    case 'pressure'
        fld = p;                            ttl = 'Pressure';
    case 'density'
        fld = rho;                          ttl = 'Density';
    case 'velocity'
        fld = hypot(u,v);                   ttl = 'Velocity magnitude';
    case 'temperature'
        fld = g*p./rho;                     ttl = 'Temperature T/T_\infty';
    case 'vorticity'
        [gvx,~] = green_gauss(Wp(:,:,3),C);
        [~,guy] = green_gauss(Wp(:,:,2),C);
        fld = (gvx(2:Ni+1,2:Nj+1) - guy(2:Ni+1,2:Nj+1))*Lref/U_inf;
        ttl = 'Vorticity \omega L_{ref}/U_\infty';
    otherwise
        error('Unknown plotField: %s',name);
end
end
