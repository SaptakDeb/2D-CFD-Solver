function ns_solver(P)
%NS_SOLVER  2D compressible Navier-Stokes / Euler, cell-centred finite volume,
%           unstructured mixed triangle / quad meshes.
%
%  Mesh        : Gmsh .msh (2.2 / 4.x ASCII), mesh2d.py text files, or built-in O-grids
%  Convective  : Rusanov or HLLC, linear reconstruction (least-squares gradients) with
%                Venkatakrishnan or Barth-Jespersen limiter
%  Viscous     : least-squares cell gradients, corrected face gradients (Re = Inf -> Euler)
%  Time        : SSP-RK2 / RK3, local (steady) or global (transient) time step
%  CFL         : adaptive (off / smart / ser) with roll-back of rejected steps
%  Wall        : slip (Euler), or no-slip adiabatic / isothermal
%  Far field   : characteristic (Riemann invariant) boundary condition
%
%  Non-dimensionalisation: rho_inf = 1, a_inf = 1, T_inf = 1, p_inf = 1/gamma, U_inf = Mach.
%  The angle of attack rotates the freestream direction, the mesh stays fixed.
%  Run through ns_input.m.

P = set_defaults(P);

gamma = P.gamma;
M_inf = P.Mach;
al    = P.aoa_deg*pi/180;
rho_inf = 1;  a_inf = 1;  p_inf = 1/gamma;
u_inf = M_inf*cos(al);
v_inf = M_inf*sin(al);
q_inf = 0.5*rho_inf*M_inf^2;
transient = strcmpi(P.mode,'transient');

%% Mesh
fprintf('\nReading mesh ...\n');
tic;
Mraw = read_mesh(P);
F = build_fvmesh(Mraw,P);
Lref = F.Lref;
fprintf('  %d cells (%d triangles, %d quads), %d edges, %d wall edges, %d far-field edges  (%.1f s)\n', ...
        F.nC,sum(F.isTri),sum(~F.isTri),F.nE,numel(F.wallE),numel(F.bFar),toc);
fprintf('  reference length %.4g, moment reference point (%.4g, %.4g), min cell area %.3e\n', ...
        Lref,F.momentRef(1),F.momentRef(2),min(F.Om));

visc = isfinite(P.Re);
wallType = lower(P.wall);
if ~visc
    wallType = 'slip';
end

C = F;
C.gamma = gamma;   C.Pr = P.Pr;
C.visc  = visc;    C.wallType = wallType;
C.mu0   = M_inf*Lref/P.Re;
C.suth  = strcmpi(P.viscosity,'sutherland');
C.Ssuth = 110.4/P.Tinf_dim;
C.cpnd  = 1/(gamma-1);
C.Tw    = P.Tw_ratio;
C.flux  = lower(P.flux);
C.limiter = lower(P.limiter);
C.venK  = P.venK;
C.rho_inf = rho_inf;  C.p_inf = p_inf;  C.a_inf = a_inf;
C.u_inf = u_inf;      C.v_inf = v_inf;
C.rho_floor = 0.01*rho_inf;
C.p_floor   = 0.01*p_inf;
C.bnx = F.nx(F.bE);
C.bny = F.ny(F.bE);

nC = F.nC;
nStages = 2 + strcmpi(P.timeScheme,'rk3');
tEnd = P.tEnd_conv*Lref/M_inf;

%% Initial condition
Uc = zeros(nC,4);
Uc(:,1) = rho_inf;
Uc(:,2) = rho_inf*u_inf;
Uc(:,3) = rho_inf*v_inf;
Uc(:,4) = p_inf/(gamma-1) + 0.5*rho_inf*M_inf^2;

if P.perturbAmp > 0 && transient
    x0 = F.xLE(1) + 1.5*Lref;
    vp = P.perturbAmp*M_inf*exp(-((F.xc(:,1)-x0).^2 + (F.xc(:,2)-F.xLE(2)-0.25*Lref).^2)/(0.5*Lref)^2);
    Uc(:,3) = Uc(:,3) + rho_inf*vp;
end

%% Solver start
fprintf('\n');
fprintf('2D Navier-Stokes Solver (unstructured, cell-centred FVM)\n');
fprintf('Mesh              : %s\n',ternary(strcmpi(P.meshSource,'file'),P.meshFile,['builtin ' P.geometry]));
fprintf('Flow              : Mach %.3f, AoA %.2f deg, %s\n',M_inf,P.aoa_deg,ternary(visc,sprintf('Re %.4g (%s viscosity, %s wall)',P.Re,P.viscosity,wallType),'Euler (slip wall)'));
fprintf('Mode              : %s\n',P.mode);
fprintf('Convective flux   : %s, order %d, limiter %s\n',C.flux,P.order,C.limiter);
fprintf('Time scheme       : SSP-%s, CFL control: %s (start %.2f, range %.2f - %.2f)\n', ...
        upper(P.timeScheme),P.cflAdapt,P.CFL,P.cflMin,P.cflMax);
fprintf('Moment            : about (%.4g, %.4g), coefficients based on Lref = %.4g\n\n', ...
        F.momentRef(1),F.momentRef(2),Lref);

maxIter = P.maxIter;
resHist = nan(maxIter,1);  ClH = nan(maxIter,1);  CdH = nan(maxIter,1);
CmH = nan(maxIter,1);      cflH = nan(maxIter,1); tH = nan(maxIter,1);
nRej = 0;
nRejRun = 0;

cfl = P.CFL;
adaptOn = ~strcmpi(P.cflAdapt,'off');
goodCount = 0;
haveStep = false;
Uprev = Uc;
resPrev = Inf;
resMax = 0;
res1 = [];
iter = 0;
t = 0;
hf = [];
if transient
    hf = figure('Name','Transient flow');
end
tic;
status = 'max iterations reached';

while iter < maxIter

    ord = ternary(iter < P.nFirst,1,P.order);

    % ---- state, residual of the current state ---------------------------------
    if ~all(isfinite(Uc(:)))
        if haveStep && adaptOn && nRejRun < 25
            Uc = Uprev;  cfl = max(P.cflMin,cfl*P.cflDown);  nRej = nRej + 1;
            nRejRun = nRejRun + 1;
            goodCount = 0;  haveStep = false;
            continue;
        end
        status = 'DIVERGED';
        break;
    end
    [W,nfl] = prim_from_cons(Uc,C);
    [Res,info] = residual(W,C,ord);
    res = sqrt(mean(Res(:,1).^2));

    guardOn = adaptOn && haveStep && cfl > P.cflMin*1.0001;
    jump = ~transient && iter > P.cflGuardIter && res > P.cflJump*resPrev;
    if guardOn && nRejRun < 25 && (nfl > 0 || ~isfinite(res) || jump)
        Uc = Uprev;
        cfl = max(P.cflMin,cfl*P.cflDown);
        nRej = nRej + 1;
        nRejRun = nRejRun + 1;
        goodCount = 0;
        haveStep = false;
        continue;
    end

    % ---- accepted state: log -----------------------------------------------------
    iter = iter + 1;
    resMax = max(resMax,res);
    if isempty(res1)
        res1 = res;
    elseif iter <= 5
        res1 = max(res1,res);
    end
    nRejRun = 0;
    resHist(iter) = res/resMax;
    [Cl,Cd,Cm,Cdp,Cdv] = coefficients(W,info,C,al,q_inf,Lref,F.momentRef);
    ClH(iter) = Cl;  CdH(iter) = Cd;  CmH(iter) = Cm;
    cflH(iter) = cfl;  tH(iter) = t;

    if mod(iter,P.printEvery) == 0 || iter == 1
        fprintf('Iter %6d  res %.3e  CFL %5.2f  Cl %9.5f  Cd %9.5f (p %.5f v %.5f)  Cm %9.5f  rej %d  t* %.3g  [%.0f s]\n', ...
                iter,resHist(iter),cfl,Cl,Cd,Cdp,Cdv,Cm,nRej,t*M_inf/Lref,toc);
    end

    % ---- CFL update ----------------------------------------------------------------
    resPrev = res;
    if adaptOn
        goodCount = goodCount + 1;
        switch lower(P.cflAdapt)
            case 'smart'
                if goodCount >= P.cflGrowEvery
                    cfl = min(P.cflMax,cfl*P.cflUp);
                    goodCount = 0;
                end
            case 'ser'
                target = P.CFL*(res1/max(res,eps))^P.cflSerExp;
                cfl = min([P.cflMax, cfl*P.cflUp^2, max(P.cflMin,target)]);
        end
    else
        cfl = P.CFL*min(1,0.2 + 0.8*iter/P.rampIter);
    end

    % ---- stopping tests ------------------------------------------------------------
    if ~transient
        if resHist(iter) < P.tol
            status = 'converged (residual)';
            break;
        end
        w = P.coefWindow;
        if iter > w + 100 && mod(iter,50) == 0
            ii = iter-w+1:iter;
            dCl = max(ClH(ii)) - min(ClH(ii));
            dCd = max(CdH(ii)) - min(CdH(ii));
            dCm = max(CmH(ii)) - min(CmH(ii));
            if dCl < P.coefTol*max(abs(mean(ClH(ii))),0.01) && ...
               dCd < P.coefTol*max(abs(mean(CdH(ii))),0.01) && ...
               dCm < P.coefTol*max(abs(mean(CmH(ii))),0.01)
                status = sprintf('coefficients stable over %d iterations',w);
                break;
            end
        end
    end
    if transient && t >= tEnd
        status = 'end time reached';
        break;
    end

    % ---- time step ---------------------------------------------------------------------
    Uprev = Uc;
    Uc0 = Uc;
    dtl = local_dt(W,C,cfl);
    if transient
        dtl = min(dtl)*ones(nC,1);
        t = t + dtl(1);
    end

    Uc = Uc0 + dtl.*Res;
    for stage = 2:nStages
        W2 = prim_from_cons(Uc,C);
        R2 = residual(W2,C,ord);
        if stage == 2
            if nStages == 3
                Uc = 0.75*Uc0 + 0.25*(Uc + dtl.*R2);
            else
                Uc = 0.5*Uc0 + 0.5*(Uc + dtl.*R2);
            end
        else
            Uc = Uc0/3 + (2/3)*(Uc + dtl.*R2);
        end
    end
    haveStep = true;

    % ---- transient plots -------------------------------------------------------------------
    if transient && mod(iter,P.plotEvery) == 0
        try
            [fld,ttl] = field_for_plot(W,info,C,P.plotField,Lref,M_inf);
            figure(hf);  clf;
            draw_field(F,fld,ttl,Lref);
            title(sprintf('%s   t* = %.2f   iter = %d',ttl,t*M_inf/Lref,iter));
            drawnow;
            if P.saveFrames
                print(hf,sprintf('frame_%05d.png',iter),'-dpng','-r120');
            end
        catch err
            fprintf('  (plot skipped: %s)\n',err.message);
        end
    end
end

resHist = resHist(1:iter);  ClH = ClH(1:iter);  CdH = CdH(1:iter);
CmH = CmH(1:iter);  cflH = cflH(1:iter);  tH = tH(1:iter);

fprintf('\nStatus: %s after %d iterations (%d rejected steps, final CFL %.2f, %.0f s)\n',status,iter,nRej,cfl,toc);
if strcmp(status,'DIVERGED')
    return;
end

%% Final state and results
[W,~] = prim_from_cons(Uc,C);
[Res,info] = residual(W,C,P.order);    %#ok<ASGLU>
[Cl,Cd,Cm,Cdp,Cdv] = coefficients(W,info,C,al,q_inf,Lref,F.momentRef);

% wall distributions
we  = F.wallE;
nxw = F.nx(we);  nyw = F.ny(we);
xw  = F.xm(we,:);
pw  = info.pw;
Cp_w = (pw - p_inf)/q_inf;
Cf_w = zeros(size(we));
if visc
    tv = info.Hv(we,2:3);                               % (tau.n) per unit area
    tn = tv(:,1).*nxw + tv(:,2).*nyw;
    tt = [tv(:,1) - tn.*nxw, tv(:,2) - tn.*nyw];
    sx = -nyw;  sy = nxw;                               % tangent, flipped to point to +x
    flip = sx < 0;
    sx(flip) = -sx(flip);  sy(flip) = -sy(flip);
    Cf_w = -(tt(:,1).*sx + tt(:,2).*sy)/q_inf;
end

% separation point on the upper surface (first sign change of Cf after the leading edge)
xsep = NaN;
if visc
    up = find(xw(:,2) > F.xLE(2) & xw(:,1) > F.xLE(1) + 0.1*Lref);
    [xs,is] = sort(xw(up,1));
    cf = Cf_w(up(is));
    k = find(cf(1:end-1) > 0 & cf(2:end) <= 0,1);
    if ~isempty(k)
        xsep = (xs(k) - cf(k)*(xs(k+1)-xs(k))/(cf(k+1)-cf(k)) - F.xLE(1))/Lref;
    end
end

fprintf('\n=== Results ===\n');
fprintf('Cl  = %.5f\n',Cl);
fprintf('Cd  = %.5f   (pressure %.5f, friction %.5f)\n',Cd,Cdp,Cdv);
fprintf('Cm  = %.5f   (about x = %.4g, positive nose-up)\n',Cm,F.momentRef(1));
if visc
    fprintf('x_sep/c (upper surface) = %.4f\n',xsep);
end
if ~transient && iter > 200
    w = min(iter-1,500);
    ii = iter-w+1:iter;
    fprintf('Variation over the last %d iterations: Cl %.2e, Cd %.2e, Cm %.2e\n',w, ...
            max(ClH(ii))-min(ClH(ii)),max(CdH(ii))-min(CdH(ii)),max(CmH(ii))-min(CmH(ii)));
end

wa = min(iter-1,P.avgWindow);
ia = iter-wa+1:iter;
ClM = mean(ClH(ia));  CdM = mean(CdH(ia));  CmM = mean(CmH(ia));
ClA = 0.5*(max(ClH(ia)) - min(ClH(ia)));
CdA = 0.5*(max(CdH(ia)) - min(CdH(ia)));
CmA = 0.5*(max(CmH(ia)) - min(CmH(ia)));
fprintf('Mean over the last %d iterations (+- half peak-to-peak):\n',wa);
fprintf('  Cl = %.5f +- %.5f,  Cd = %.5f +- %.5f,  Cm = %.5f +- %.5f\n',ClM,ClA,CdM,CdA,CmM,CmA);

results = struct('Cl',Cl,'Cd',Cd,'Cdp',Cdp,'Cdv',Cdv,'Cm',Cm,'xsep',xsep,'iter',iter,'status',status, ...
                 'ClMean',ClM,'CdMean',CdM,'CmMean',CmM,'ClAmp',ClA,'CdAmp',CdA,'CmAmp',CmA, ...
                 'ClHist',ClH,'CdHist',CdH,'CmHist',CmH,'resHist',resHist, ...
                 'xw',xw(:,1),'yw',xw(:,2),'Cp',Cp_w,'Cf',Cf_w);          %#ok<NASGU>
assignin('base','ns_results',results);

%% Plots
if P.plots
    try
        vort = info.gx(1:nC,3) - info.gy(1:nC,2);
        Mach = hypot(W(1:nC,2),W(1:nC,3))./sqrt(gamma*W(1:nC,4)./W(1:nC,1));
        Cp   = (W(1:nC,4) - p_inf)/q_inf;
        wx = F.xLE(1) + [-0.3 1.6]*Lref;
        wy = F.xLE(2) + [-0.5 0.5]*Lref;

        figure;  draw_field(F,Mach,'Mach number',Lref);  title('Mach number');
        figure;  draw_field(F,Cp,'C_p',Lref);            title('Pressure coefficient');
        figure;  draw_field(F,vort*Lref/M_inf,'vorticity',Lref);
        vm = 0.3*max(abs(vort*Lref/M_inf));  caxis([-vm vm]);  title('Vorticity \omega L_{ref}/U_\infty');
        figure;  draw_mesh(F,wx,wy);                     title('Mesh near the body');

        figure;
        subplot(2,1,1);
        plot_surface(xw,Cp_w,F.xLE(2));  set(gca,'YDir','reverse');  grid on;
        xlabel('x'); ylabel('C_p'); title('Wall pressure coefficient');
        subplot(2,1,2);
        plot_surface(xw,Cf_w,F.xLE(2));  grid on;
        xlabel('x'); ylabel('C_f'); title('Skin friction');

        figure;
        xh = ternary(transient,tH*M_inf/Lref,(1:iter)');
        xlab = ternary(transient,'t^* = t U_\infty / L_{ref}','iteration');
        subplot(4,1,1); semilogy(xh,resHist,'LineWidth',1.5); grid on;
        ylabel('density residual'); title('Convergence history');
        subplot(4,1,2); plot(xh,ClH,xh,CdH,'LineWidth',1.5); grid on; legend('C_l','C_d');
        subplot(4,1,3); plot(xh,CmH,'LineWidth',1.5); grid on; ylabel('C_m');
        subplot(4,1,4); plot(xh,cflH,'LineWidth',1.5); grid on; ylabel('CFL'); xlabel(xlab);
    catch err
        fprintf('(plots skipped: %s)\n',err.message);
    end
end

end


%% ========================================================================
%  Subfunctions
%  ========================================================================

function out = ternary(c,a,b)
if c
    out = a;
else
    out = b;
end
end


function P = set_defaults(P)
D.meshSource = 'builtin';   D.meshFile = '';
D.geometry   = 'cylinder';
D.Rc = 0.5;  D.Lb = 1.0;  D.Hb = 0.7;  D.naca_t = 0.12;
D.Ni = 160;  D.Nj = 90;  D.h1_Lref = 0.015;  D.Rfar_Lref = 20;  D.blendPower = 1.5;
D.wallTags = [];  D.farTags = [];  D.Lref = [];  D.momentRef = [];
D.aoa_deg = 0;
D.Mach = 0.3;  D.Re = 40;  D.Pr = 0.72;  D.gamma = 1.4;
D.viscosity = 'constant';  D.Tinf_dim = 288.15;
D.wall = 'adiabatic';  D.Tw_ratio = 1.0;
D.flux = 'hllc';  D.order = 2;  D.limiter = 'venkat';  D.venK = 3;
D.timeScheme = 'rk3';  D.nFirst = 200;  D.rampIter = 300;
D.CFL = 1.0;  D.cflAdapt = 'smart';
D.cflMin = 0.1;  D.cflMax = 3.0;  D.cflUp = 1.05;  D.cflDown = 0.5;
D.cflGrowEvery = 10;  D.cflJump = 5;  D.cflGuardIter = 50;  D.cflSerExp = 1.0;
D.mode = 'steady';  D.maxIter = 30000;  D.tol = 1e-8;
D.coefWindow = 1000;  D.coefTol = 2e-3;  D.avgWindow = 3000;
D.tEnd_conv = 100;  D.plotEvery = 500;  D.plotField = 'vorticity';
D.saveFrames = false;  D.perturbAmp = 1e-2;
D.printEvery = 200;  D.plots = true;
fn = fieldnames(D);
for k = 1:numel(fn)
    if ~isfield(P,fn{k})
        P.(fn{k}) = D.(fn{k});
    end
end
end


function [W,nfl] = prim_from_cons(Uc,C)
% primitive variables [rho u v p] on cells + ghost cells; nfl = number of floored cells
g = C.gamma;
rho0 = Uc(:,1);
rho = max(rho0,C.rho_floor);
u = Uc(:,2)./rho;
v = Uc(:,3)./rho;
p0 = (g-1)*(Uc(:,4) - 0.5*rho.*(u.^2 + v.^2));
p = max(p0,C.p_floor);
nfl = sum(rho0 <= C.rho_floor | p0 <= C.p_floor);
W = [rho u v p; zeros(C.nB,4)];
W = apply_bc(W,C);
end


function W = apply_bc(W,C)
nC = C.nC;  g = C.gamma;

% ---- wall ghost cells -------------------------------------------------------------------
k  = C.bWall;
Lc = C.bL(k);
nx = C.bnx(k);  ny = C.bny(k);
r = W(Lc,1);  u = W(Lc,2);  v = W(Lc,3);  p = W(Lc,4);
switch C.wallType
    case 'slip'
        vn = u.*nx + v.*ny;
        ug = u - 2*vn.*nx;   vg = v - 2*vn.*ny;   rg = r;
    case 'isothermal'
        T  = g*p./r;
        Tg = max(2*C.Tw - T,0.2*C.Tw);
        ug = -u;   vg = -v;   rg = g*p./Tg;
    otherwise                                      % adiabatic no-slip
        ug = -u;   vg = -v;   rg = r;
end
W(nC+k,:) = [rg ug vg p];

% ---- far-field ghost cells (characteristic) ----------------------------------------------------
k  = C.bFar;
Lc = C.bL(k);
nx = C.bnx(k);  ny = C.bny(k);
rI = W(Lc,1);  uI = W(Lc,2);  vI = W(Lc,3);  pI = W(Lc,4);
aI  = sqrt(g*pI./rI);
unI = uI.*nx + vI.*ny;
unF = C.u_inf*nx + C.v_inf*ny;
Rp  = unI + 2*aI/(g-1);
Rm  = unF - 2*C.a_inf/(g-1);
unB = 0.5*(Rp + Rm);
aB  = max(0.25*(g-1)*(Rp - Rm),0.1*C.a_inf);
isOut = double(unB > 0);
sI  = pI./rI.^g;
sF  = C.p_inf/C.rho_inf^g;
sB  = isOut.*sI + (1-isOut)*sF;
utx = isOut.*(uI - unI.*nx) + (1-isOut).*(C.u_inf - unF.*nx);
uty = isOut.*(vI - unI.*ny) + (1-isOut).*(C.v_inf - unF.*ny);
rB  = (aB.^2./(g*sB)).^(1/(g-1));
pB  = rB.*aB.^2/g;
W(nC+k,:) = [rB, unB.*nx + utx, unB.*ny + uty, pB];
end


function mu = viscosity(T,C)
if C.suth
    mu = C.mu0*T.^1.5.*(1 + C.Ssuth)./(T + C.Ssuth);
else
    mu = C.mu0*ones(size(T));
end
end


function dtl = local_dt(W,C,cfl)
g = C.gamma;
nC = C.nC;
rho = W(:,1);  u = W(:,2);  v = W(:,3);  p = W(:,4);
a  = sqrt(g*p./rho);
Lc_e = (abs(0.5*(u(C.L)+u(C.R)).*C.nx + 0.5*(v(C.L)+v(C.R)).*C.ny) ...
        + 0.5*(a(C.L)+a(C.R))).*C.A;
Lc = C.ML*Lc_e + C.MR*Lc_e;
Lv = zeros(nC,1);
if C.visc
    mu = viscosity(g*p./rho,C);
    nu = max(4/3,g/C.Pr)*mu./rho;
    nu_e = 0.5*(nu(C.L) + nu(C.R));
    Lv_e = nu_e.*C.A.^2;
    Lv = (C.ML*Lv_e + C.MR*Lv_e)./C.Om;
end
dtl = cfl*C.Om./(Lc + 4*Lv);
end


function F = conv_flux(WL,WR,nx,ny,g,flux)
rL = WL(:,1);  uL = WL(:,2);  vL = WL(:,3);  pL = WL(:,4);
rR = WR(:,1);  uR = WR(:,2);  vR = WR(:,3);  pR = WR(:,4);

unL =  uL.*nx + vL.*ny;   utL = -uL.*ny + vL.*nx;
unR =  uR.*nx + vR.*ny;   utR = -uR.*ny + vR.*nx;

EL = pL/(g-1) + 0.5*rL.*(unL.^2 + utL.^2);
ER = pR/(g-1) + 0.5*rR.*(unR.^2 + utR.^2);

UL = [rL, rL.*unL, rL.*utL, EL];
UR = [rR, rR.*unR, rR.*utR, ER];
FL = [rL.*unL, rL.*unL.^2 + pL, rL.*unL.*utL, unL.*(EL + pL)];
FR = [rR.*unR, rR.*unR.^2 + pR, rR.*unR.*utR, unR.*(ER + pR)];

aL = sqrt(g*pL./rL);
aR = sqrt(g*pR./rR);

if strcmp(flux,'rusanov')
    smax = max(abs(unL) + aL,abs(unR) + aR);
    Fn = 0.5*(FL + FR) - 0.5*smax.*(UR - UL);
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
    UsL = [kL, kL.*Ss, kL.*utL, kL.*(EL./rL + (Ss - unL).*(Ss + pL./(rL.*(SL - unL))))];
    UsR = [kR, kR.*Ss, kR.*utR, kR.*(ER./rR + (Ss - unR).*(Ss + pR./(rR.*(SR - unR))))];
    FsL = FL + SL.*(UsL - UL);
    FsR = FR + SR.*(UsR - UR);

    m1 = double(SL >= 0);
    m2 = double(SL <  0 & Ss >= 0);
    m3 = double(Ss <  0 & SR >  0);
    m4 = double(SR <= 0);
    Fn = m1.*FL + m2.*FsL + m3.*FsR + m4.*FR;
end

F = [Fn(:,1), Fn(:,2).*nx - Fn(:,3).*ny, Fn(:,2).*ny + Fn(:,3).*nx, Fn(:,4)];
end


function [R,info] = residual(W,C,ord)
g  = C.gamma;
nC = C.nC;
L  = C.L;
Rr = C.R;

rho = W(:,1);  u = W(:,2);  v = W(:,3);  p = W(:,4);
T   = g*p./rho;
Phi = [rho u v p T];

% ---- least-squares gradients (cells), ghost cells copy the gradient of their cell ------------
dPhi = Phi(Rr,:) - Phi(L,:);
bx = C.ML*(C.wdx.*dPhi) + C.MR*(C.wdx.*dPhi);
by = C.ML*(C.wdy.*dPhi) + C.MR*(C.wdy.*dPhi);
gx = C.Mi11.*bx + C.Mi12.*by;
gy = C.Mi12.*bx + C.Mi22.*by;
gx = [gx; gx(C.bL,:)];
gy = [gy; gy(C.bL,:)];

% ---- limited linear reconstruction -----------------------------------------------------------------
psi = zeros(nC + C.nB,4);
if ord == 2
    phiC = Phi(1:nC,1:4);
    pmax = phiC;  pmin = phiC;
    for j = 1:4
        nbv  = Phi(C.NB(:,j),1:4);
        pmax = max(pmax,nbv);
        pmin = min(pmin,nbv);
    end
    ps = ones(nC,4);
    eps2 = (C.venK*C.h).^3;
    for j = 1:4
        d2 = gx(1:nC,1:4).*C.dxf(:,j) + gy(1:nC,1:4).*C.dyf(:,j);
        d1 = (pmax - phiC).*(d2 > 0) + (pmin - phiC).*(d2 <= 0);
        if strcmp(C.limiter,'barth')
            r = d1./d2;
            r(d2 == 0) = 1;
            pj = min(1,r);
        else
            pj = ((d1.^2 + eps2) + 2*d2.*d1)./(d1.^2 + 2*d2.^2 + d1.*d2 + eps2);
        end
        ps = min(ps,pj);
    end
    psi(1:nC,:) = max(min(ps,1),0);
end

WL = Phi(L,1:4)  + psi(L,:).*(gx(L,1:4).*C.dLx + gy(L,1:4).*C.dLy);
WR = Phi(Rr,1:4) + psi(Rr,:).*(gx(Rr,1:4).*C.dRx + gy(Rr,1:4).*C.dRy);
bad = WL(:,1) <= 0 | WL(:,4) <= 0 | WR(:,1) <= 0 | WR(:,4) <= 0;
if any(bad)
    WL(bad,:) = Phi(L(bad),1:4);
    WR(bad,:) = Phi(Rr(bad),1:4);
end

% ---- convective flux -----------------------------------------------------------------------------------------
Fc = conv_flux(WL,WR,C.nx,C.ny,g,C.flux);
we = C.wallE;
pw = WL(we,4);
Fc(we,:) = [zeros(numel(we),1), pw.*C.nx(we), pw.*C.ny(we), zeros(numel(we),1)];

% ---- viscous flux -------------------------------------------------------------------------------------------------
Hv = zeros(C.nE,4);
if C.visc
    mu = viscosity(T,C);
    kk = mu*C.cpnd/C.Pr;
    ae = 0.5*(Phi(L,:) + Phi(Rr,:));
    mue = 0.5*(mu(L) + mu(Rr));
    ke  = 0.5*(kk(L) + kk(Rr));
    corr = (Phi(Rr,:) - Phi(L,:))./C.ds ...
           - (0.5*(gx(L,:)+gx(Rr,:)).*C.ex + 0.5*(gy(L,:)+gy(Rr,:)).*C.ey);
    fx = 0.5*(gx(L,:) + gx(Rr,:)) + corr.*C.ex;
    fy = 0.5*(gy(L,:) + gy(Rr,:)) + corr.*C.ey;
    ux = fx(:,2);  uy = fy(:,2);
    vx = fx(:,3);  vy = fy(:,3);
    Tx = fx(:,5);  Ty = fy(:,5);
    div = ux + vy;
    txx = mue.*(2*ux - (2/3)*div);
    tyy = mue.*(2*vy - (2/3)*div);
    txy = mue.*(uy + vx);
    Hx = txx.*C.nx + txy.*C.ny;
    Hy = txy.*C.nx + tyy.*C.ny;
    He = (ae(:,2).*txx + ae(:,3).*txy + ke.*Tx).*C.nx + (ae(:,2).*txy + ae(:,3).*tyy + ke.*Ty).*C.ny;
    Hv = [zeros(C.nE,1), Hx, Hy, He];
end

H = (Fc - Hv).*C.A;
R = -(C.ML*H - C.MR*H)./C.Om;

info.pw = pw;
info.Hv = Hv;
info.gx = gx;
info.gy = gy;
end


function [Cl,Cd,Cm,Cdp,Cdv] = coefficients(W,info,C,al,q_inf,Lref,xref) %#ok<INUSL>
we = C.wallE;
nx = C.nx(we);  ny = C.ny(we);  A = C.A(we);
% pressure pushes the body along the fluid-outward normal, viscous traction acts opposite to tau.n
Fpx =  info.pw.*nx.*A;   Fpy =  info.pw.*ny.*A;
Fvx = -info.Hv(we,2).*A; Fvy = -info.Hv(we,3).*A;
xm = C.xm(we,:);
rx = xm(:,1) - xref(1);
ry = xm(:,2) - xref(2);
Fx = Fpx + Fvx;  Fy = Fpy + Fvy;
ca = cos(al);  sa = sin(al);
Cd  = (sum(Fx)*ca + sum(Fy)*sa)/(q_inf*Lref);
Cl  = (-sum(Fx)*sa + sum(Fy)*ca)/(q_inf*Lref);
Cdp = (sum(Fpx)*ca + sum(Fpy)*sa)/(q_inf*Lref);
Cdv = (sum(Fvx)*ca + sum(Fvy)*sa)/(q_inf*Lref);
Cm  = sum(rx.*Fy - ry.*Fx)/(q_inf*Lref^2);
end


function [fld,ttl] = field_for_plot(W,info,C,name,Lref,M_inf)
nC = C.nC;  g = C.gamma;
rho = W(1:nC,1);  u = W(1:nC,2);  v = W(1:nC,3);  p = W(1:nC,4);
switch lower(name)
    case 'mach',        fld = hypot(u,v)./sqrt(g*p./rho);   ttl = 'Mach number';
    case 'pressure',    fld = p;                            ttl = 'Pressure';
    case 'density',     fld = rho;                          ttl = 'Density';
    case 'velocity',    fld = hypot(u,v);                   ttl = 'Velocity magnitude';
    case 'temperature', fld = g*p./rho;                     ttl = 'Temperature';
    case 'vorticity',   fld = (info.gx(1:nC,3) - info.gy(1:nC,2))*Lref/M_inf;  ttl = 'Vorticity';
    otherwise,          error('unknown plotField %s',name);
end
end


function draw_field(F,fld,ttl,Lref)
hold on;
t = F.isTri;
tri = F.cells(t,1:3);
qd  = F.cells(~t,:);
if ~isempty(tri)
    patch('Faces',tri,'Vertices',F.nodes,'FaceVertexCData',fld(t),'FaceColor','flat','EdgeColor','none');
end
if ~isempty(qd)
    patch('Faces',qd,'Vertices',F.nodes,'FaceVertexCData',fld(~t),'FaceColor','flat','EdgeColor','none');
end
we = F.wallE;
axis equal;
xlim(F.xLE(1) + [-0.3 1.6]*Lref);
ylim(F.xLE(2) + [-0.5 0.5]*Lref);
colorbar;
xlabel('x');  ylabel('y');
if strcmpi(ttl,'vorticity') || strcmpi(ttl,'Vorticity')
    vm = 0.3*max(abs(fld));
    caxis([-vm vm]);
end
end


function draw_mesh(F,wx,wy)
hold on;
t = F.isTri;
if any(t)
    patch('Faces',F.cells(t,1:3),'Vertices',F.nodes,'FaceColor','none','EdgeColor','k','LineWidth',0.3);
end
if any(~t)
    patch('Faces',F.cells(~t,:),'Vertices',F.nodes,'FaceColor','none','EdgeColor','k','LineWidth',0.3);
end
axis equal;
xlim(wx);  ylim(wy);
end


function plot_surface(xw,val,yref)
up = xw(:,2) >= yref;
[xa,ia] = sort(xw(up,1));
va = val(up);
[xb,ib] = sort(xw(~up,1));
vb = val(~up);
plot(xa,va(ia),'b-',xb,vb(ib),'r-','LineWidth',1.3);
legend('upper','lower');
end
