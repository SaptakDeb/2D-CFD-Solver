function sol = euler_solver(cfg, m)
% EULER_SOLVER  2D compressible Euler equations, cell-centred finite volume
% on a body-fitted quad mesh (any topology produced by the mesh_*.m files).
%
%   - face-based fluxes: Rusanov / HLL / HLLC
%   - MUSCL reconstruction from least-squares gradients, with the TVD
%     limiters minmod / vanleer / mc / vanalbada / superbee
%   - SSP Runge-Kutta 2 or 3, steady (local dt) or transient (global dt)
%   - slip wall by pressure flux on the exact surface, characteristic far field

fl = cfg.flow;  nm = cfg.num;
gam = fl.gamma;  R = fl.R;
aoa = fl.aoa_deg*pi/180;
rho_inf = fl.p_inf/(R*fl.T_inf);
a_inf   = sqrt(gam*R*fl.T_inf);
U_inf   = fl.M_inf*a_inf;
u_inf   = U_inf*cos(aoa);   v_inf = U_inf*sin(aoa);   % flow in the body frame
q_inf   = 0.5*rho_inf*U_inf^2;
p_inf   = fl.p_inf;
Lref    = m.Lref;
transient = strcmpi(nm.mode,'transient');
tEnd = nm.tEnd_conv*Lref/U_inf;

Nc = numel(m.vol);  nIF = numel(m.iL);  nBF = numel(m.bL);

%% Operators and constants used by the residual
S = struct();
S.iL = m.iL;  S.iR = m.iR;  S.inx = m.inx;  S.iny = m.iny;  S.ilen = m.ilen;
S.vol = m.vol;  S.nb0 = m.nb0;  S.cx = m.cx;  S.cy = m.cy;
S.bL = m.bL;  S.bnx = m.bnx;  S.bny = m.bny;  S.blen = m.blen;
S.Sint = sparse([m.iL; m.iR],[(1:nIF)'; (1:nIF)'],[-ones(nIF,1); ones(nIF,1)],Nc,nIF);
S.Sbnd = sparse(m.bL,(1:nBF)',-ones(nBF,1),Nc,nBF);
S.AL = sparse(m.iL,(1:nIF)',ones(nIF,1),Nc,nIF);
S.AR = sparse(m.iR,(1:nIF)',ones(nIF,1),Nc,nIF);
S.AB = sparse(m.bL,(1:nBF)',ones(nBF,1),Nc,nBF);

dxLR = m.xc(m.iR) - m.xc(m.iL);   dyLR = m.yc(m.iR) - m.yc(m.iL);
dd   = dxLR.^2 + dyLR.^2;
S.dxLR = dxLR;  S.dyLR = dyLR;
S.lamL = min(max(((m.ixm - m.xc(m.iL)).*dxLR + (m.iym - m.yc(m.iL)).*dyLR)./dd, 0.1), 0.9);
S.lamR = min(max(((m.xc(m.iR) - m.ixm).*dxLR + (m.yc(m.iR) - m.iym).*dyLR)./dd, 0.1), 0.9);

S.iW = find(m.btype == 1);   S.iF = find(m.btype == 2);
S.wc  = m.bL(S.iW);  S.wnx = m.bnx(S.iW);  S.wny = m.bny(S.iW);  S.wlen = m.blen(S.iW);
S.wdx = m.bxm(S.iW) - m.xc(S.wc);   S.wdy = m.bym(S.iW) - m.yc(S.wc);
S.fc  = m.bL(S.iF);  S.fnx = m.bnx(S.iF);  S.fny = m.bny(S.iF);
S.nBF = nBF;

P = struct('gam',gam,'rho_floor',0.01*rho_inf,'p_floor',0.01*p_inf, ...
           'limiter',lower(nm.limiter),'flux',lower(nm.flux),'farfield',lower(nm.farfield), ...
           'rho_inf',rho_inf,'u_inf',u_inf,'v_inf',v_inf,'p_inf',p_inf,'a_inf',a_inf);

nStages = 2 + strcmpi(nm.timeScheme,'rk3');

%% Initial condition: freestream + small wake perturbation
rho = rho_inf*ones(Nc,1);  u = u_inf*ones(Nc,1);  v = v_inf*ones(Nc,1);  p = p_inf*ones(Nc,1);
xp = m.xb + 1.5*Lref;   yp = m.yb + 0.25*Lref;
if nm.perturbAmp > 0
    v = v + nm.perturbAmp*U_inf*exp(-((m.xc-xp).^2 + (m.yc-yp).^2)/(0.5*Lref)^2);
end
Uc = [rho, rho.*u, rho.*v, p/(gam-1) + 0.5*rho.*(u.^2 + v.^2)];
[~,ip] = min((m.xc-xp).^2 + (m.yc-yp).^2);       % wake probe cell

%% Header
dtEst = min(2*nm.CFL*m.vol./ (S.AL*((a_inf+U_inf)*m.ilen) + S.AR*((a_inf+U_inf)*m.ilen) + S.AB*((a_inf+U_inf)*m.blen)));
fprintf('\n2D Euler Solver (Version 4, body-fitted)\n');
fprintf('Body              : %s (%s, Lref = %.3f)\n',m.name,m.topology,Lref);
fprintf('Mode              : %s\n',nm.mode);
fprintf('Riemann solver    : %s\n',nm.flux);
fprintf('Reconstruction    : order %d, limiter %s (first %d iterations first order)\n',nm.order,nm.limiter,nm.nFirst);
fprintf('Time scheme       : SSP-%s, CFL %.2f\n',upper(nm.timeScheme),nm.CFL);
fprintf('Far field         : %s\n',nm.farfield);
fprintf('Mesh              : %d cells, %d wall faces\n',Nc,numel(S.iW));
fprintf('Mach / AoA        : %.3f / %.2f deg\n',fl.M_inf,fl.aoa_deg);
if transient
    fprintf('End time          : t* = %.1f  (smallest stable dt = %.2e s -> about %d steps)\n', ...
            nm.tEnd_conv,dtEst,ceil(tEnd/dtEst));
end
fprintf('\n');

%% Live plot set-up
if transient
    xl = [m.xb - 3*Lref, m.xb + 12*Lref];  yl = [m.yb - 3*Lref, m.yb + 3*Lref];
    win = find(m.xc >= xl(1) & m.xc <= xl(2) & m.yc >= yl(1) & m.yc <= yl(2));
    hf = figure('Name','Transient flow');
    hp = patch('Faces',m.cells(win,:),'Vertices',[m.xn m.yn],'FaceVertexCData',zeros(numel(win),1), ...
               'FaceColor','flat','EdgeColor','none');
    hold on;  line(m.wx,m.wy,'Color','k','LineWidth',1.5);
    axis equal;  xlim(xl);  ylim(yl);  colorbar;  xlabel('x');  ylabel('y');
end

resHist = nan(nm.maxIter,1);  ClHist = resHist;  CdHist = resHist;
tHist = resHist;  vpHist = resHist;
res0 = [];  t = 0;  tic;

%% Time marching
for iter = 1:nm.maxIter

    if iter <= nm.nFirst, ord = 1; else ord = nm.order; end
    CFLn = nm.CFL*min(1, 0.2 + 0.8*iter/nm.rampIter);
    Uc0  = Uc;

    for stage = 1:nStages
        [Res,W,aux] = residual(Uc,S,P,ord);

        if stage == 1
            dtl = CFLn*aux.dt1;
            if transient, dtl = min(dtl)*ones(Nc,1); end
            dtmin = min(dtl);
            Res1 = Res;   W1 = W;   aux1 = aux;
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

    % Histories (state at the start of the step)
    tHist(iter) = t;
    if transient, t = t + dtmin; end
    res = sqrt(mean(Res1(:,1).^2));
    if isempty(res0), res0 = res; end
    resHist(iter) = res/res0;

    Fx = aux1.Fx;  Fy = aux1.Fy;
    CdHist(iter) = (Fx*cos(aoa) + Fy*sin(aoa))/(q_inf*Lref);
    ClHist(iter) = (-Fx*sin(aoa) + Fy*cos(aoa))/(q_inf*Lref);
    vpHist(iter) = W1(ip,3);

    if ~all(isfinite(Uc(:)))
        fprintf('\nSOLVER DIVERGED\nIteration : %d\nTime      : %.2f s\n',iter,toc);
        break;
    end

    if mod(iter,nm.printEvery) == 0
        fprintf('Iteration: %6d / %6d   Residual: %.3e   Cl: %8.4f   Cd: %8.4f   t*: %7.2f   Time: %.1f s\n', ...
                iter,nm.maxIter,resHist(iter),ClHist(iter),CdHist(iter),tHist(iter)*U_inf/Lref,toc);
    end

    if transient && (mod(iter,nm.plotEvery) == 0 || iter == 1)
        [Gx,Gy] = cell_gradient(W1,S);
        [fld,ttl] = field_of(nm.plotField,W1,Gx,Gy,P,Lref,U_inf);
        set(hp,'FaceVertexCData',fld(win));
        if any(strcmpi(nm.plotField,{'vorticity'}))
            vmax = max(abs(fld(win)));  caxis([-0.3*vmax 0.3*vmax]);
        elseif strcmpi(nm.plotField,'schlieren')
            colormap(gray);  caxis([0 1]);
        else
            caxis([min(fld(win)) max(fld(win))]);
        end
        title(sprintf('%s   t* = %.2f   iter = %d',ttl,t*U_inf/Lref,iter));
        drawnow;
        if nm.saveFrames
            print(hf,sprintf('frame_%05d.png',iter),'-dpng','-r120');
        end
    end

    if transient && t >= tEnd
        fprintf('\nEnd time reached\nIterations : %d\nt*         : %.2f\nTime       : %.2f s\n',iter,t*U_inf/Lref,toc);
        break;
    end
    if ~transient && resHist(iter) < nm.tol
        fprintf('\nConverged\nIterations : %d\nResidual   : %.3e\nTime       : %.2f s\n',iter,resHist(iter),toc);
        break;
    end
end

if iter == nm.maxIter && ((~transient && resHist(iter) >= nm.tol) || (transient && t < tEnd))
    fprintf('\nMaximum iterations reached\nResidual: %.3e\nTime    : %.2f s\n',resHist(iter),toc);
end

%% Final state and fields
[~,W,aux] = residual(Uc,S,P,nm.order);
[Gx,Gy] = cell_gradient(W,S);
rho = W(:,1);  u = W(:,2);  v = W(:,3);  p = W(:,4);
V    = hypot(u,v);
Mach = V./sqrt(gam*p./rho);
Cp   = (p - p_inf)/q_inf;
dS   = p./p_inf.*(rho_inf./rho).^gam - 1;
vort = (Gx(:,3) - Gy(:,2))*Lref/U_inf;

sol.W = W;  sol.Mach = Mach;  sol.Cp = Cp;  sol.dS = dS;  sol.vort = vort;  sol.V = V;
sol.resHist = resHist(1:iter);  sol.ClHist = ClHist(1:iter);  sol.CdHist = CdHist(1:iter);
sol.tHist = tHist(1:iter);      sol.vpHist = vpHist(1:iter);
sol.iter = iter;  sol.transient = transient;  sol.U_inf = U_inf;  sol.q_inf = q_inf;
sol.wallX = m.bxm(S.iW);  sol.wallY = m.bym(S.iW);
sol.wallCp = (aux.pw - p_inf)/q_inf;
sol.Cl = (-aux.Fx*sin(aoa) + aux.Fy*cos(aoa))/(q_inf*Lref);
sol.Cd = ( aux.Fx*cos(aoa) + aux.Fy*sin(aoa))/(q_inf*Lref);

dS_L2 = sqrt(mean(dS.^2));
fprintf('\nSolution complete\n');
fprintf('Body                : %s, AoA = %.2f deg\n',m.name,fl.aoa_deg);
fprintf('Freestream velocity : %.3f m/s\n',U_inf);
fprintf('Maximum velocity    : %.3f m/s\n',max(V));
fprintf('Minimum pressure    : %.2f Pa\n',min(p));
fprintf('Maximum pressure    : %.2f Pa\n',max(p));
fprintf('Maximum Mach number : %.4f\n',max(Mach));
fprintf('Lift coefficient    : %.4f  (based on Lref = %.3f)\n',sol.Cl,Lref);
fprintf('Drag coefficient    : %.4f\n',sol.Cd);
fprintf('Entropy error L2    : %.3e\n',dS_L2);
fprintf('Entropy error max   : %.3e\n',max(dS));

%% Vortex shedding frequency (transient)
sol.St = NaN;
if transient
    th = sol.tHist;
    i2 = find(th >= 0.4*th(end));
    if numel(i2) > 64
        sigs = {sol.ClHist, sol.vpHist};  names = {'Cl signal','wake probe v'};
        tu = linspace(th(i2(1)),th(end),numel(i2));
        fs = 1/(tu(2) - tu(1));  n = numel(tu);
        win_ = 0.5*(1 - cos(2*pi*(0:n-1)/(n-1)));
        for ks = 1:2
            sg = interp1(th(i2),sigs{ks}(i2),tu);
            sg = (sg(:).' - mean(sg)).*win_;
            Sf = abs(fft(sg));
            [~,im] = max(Sf(2:floor(n/2)));
            f = im*fs/n;
            fprintf('Shedding (%s): f = %.2f Hz, St = f Lref/U_inf = %.3f\n',names{ks},f,f*Lref/U_inf);
            if ks == 2, sol.St = f*Lref/U_inf; end
        end
        fprintf('Mean Cl (last 60%%) : %.4f\n',mean(sol.ClHist(i2)));
        fprintf('Mean Cd (last 60%%) : %.4f\n',mean(sol.CdHist(i2)));
        fprintf('(Spectrum resolution is limited by the record length; run longer for a sharper St.)\n');
    else
        fprintf('Run too short for a shedding frequency estimate.\n');
    end
end
end

%% ======================================================================
function [Res,W,aux] = residual(Uc,S,P,ord)
% Flux balance dU/dt for all cells, plus wall forces and local time step.
gam = P.gam;
rho = max(Uc(:,1),P.rho_floor);
u = Uc(:,2)./rho;   v = Uc(:,3)./rho;
p = max((gam-1)*(Uc(:,4) - 0.5*rho.*(u.^2 + v.^2)),P.p_floor);
W = [rho u v p];

if ord == 2, [Gx,Gy] = cell_gradient(W,S); end

% ---- interior faces
WL = W(S.iL,:);  WR = W(S.iR,:);
if ord == 2
    dW  = WR - WL;
    dWs = dW + (dW == 0)*1e-30;
    rL = 2*(Gx(S.iL,:).*S.dxLR + Gy(S.iL,:).*S.dyLR)./dWs - 1;
    rR = 2*(Gx(S.iR,:).*S.dxLR + Gy(S.iR,:).*S.dyLR)./dWs - 1;
    rL = max(min(rL,1e3),-1e3);   rR = max(min(rR,1e3),-1e3);
    WLf = WL + S.lamL.*limiter(rL,P.limiter).*dW;
    WRf = WR - S.lamR.*limiter(rR,P.limiter).*dW;
    bad = WLf(:,1) <= 0 | WLf(:,4) <= 0 | WRf(:,1) <= 0 | WRf(:,4) <= 0;
    WLf(bad,:) = WL(bad,:);   WRf(bad,:) = WR(bad,:);
else
    WLf = WL;  WRf = WR;
end
F = numflux(WLf,WRf,S.inx,S.iny,gam,P.flux);
Res = S.Sint*(F.*S.ilen);

% ---- boundary faces
Fb = zeros(S.nBF,4);

% Wall: only a pressure flux acts. Pressure is extrapolated to the face
% and kept inside the range of the neighbouring cells.
pc = p(S.wc);
if ord == 2
    pw = pc + Gx(S.wc,4).*S.wdx + Gy(S.wc,4).*S.wdy;
    pn = [p(S.nb0(S.wc,:)), pc];
    pw = min(max(pw,min(pn,[],2)),max(pn,[],2));
else
    pw = pc;
end
Fb(S.iW,:) = [zeros(size(pw)), pw.*S.wnx, pw.*S.wny, zeros(size(pw))];

% Far field
c = S.fc;  nx = S.fnx;  ny = S.fny;
rI = rho(c);  uI = u(c);  vI = v(c);  pI = p(c);
unF = P.u_inf*nx + P.v_inf*ny;
if strcmp(P.farfield,'characteristic')
    aI  = sqrt(gam*pI./rI);
    unI = uI.*nx + vI.*ny;
    Rp  = unI + 2*aI/(gam-1);                  % outgoing invariant
    Rm  = unF - 2*P.a_inf/(gam-1);             % incoming invariant
    unB = 0.5*(Rp + Rm);
    aB  = max(0.25*(gam-1)*(Rp - Rm),0.1*P.a_inf);
    out = double(unB > 0);
    sB  = out.*(pI./rI.^gam) + (1-out)*(P.p_inf/P.rho_inf^gam);
    utx = out.*(uI - unI.*nx) + (1-out).*(P.u_inf - unF.*nx);
    uty = out.*(vI - unI.*ny) + (1-out).*(P.v_inf - unF.*ny);
    rB  = (aB.^2./(gam*sB)).^(1/(gam-1));
    pB  = rB.*aB.^2/gam;
    uB  = unB.*nx + utx;   vB = unB.*ny + uty;
else   % 'freestream': Dirichlet on inflow, extrapolation with p_inf on outflow
    out = unF > 0;
    rB = P.rho_inf*ones(size(c));  uB = P.u_inf*ones(size(c));
    vB = P.v_inf*ones(size(c));    pB = P.p_inf*ones(size(c));
    rB(out) = rI(out);  uB(out) = uI(out);  vB(out) = vI(out);
end
Fb(S.iF,:) = physflux(rB,uB,vB,pB,nx,ny,gam);

Res = (Res + S.Sbnd*(Fb.*S.blen))./S.vol;

% ---- wall force (pressure acting on the body; wall normal points into the body)
aux.Fx = sum((pw - P.p_inf).*S.wnx.*S.wlen);
aux.Fy = sum((pw - P.p_inf).*S.wny.*S.wlen);
aux.pw = pw;

% ---- local time step for CFL = 1
a  = sqrt(gam*p./rho);
lL = (abs(u(S.iL).*S.inx + v(S.iL).*S.iny) + a(S.iL)).*S.ilen;
lR = (abs(u(S.iR).*S.inx + v(S.iR).*S.iny) + a(S.iR)).*S.ilen;
lB = (abs(u(S.bL).*S.bnx + v(S.bL).*S.bny) + a(S.bL)).*S.blen;
aux.dt1 = 2*S.vol./(S.AL*lL + S.AR*lR + S.AB*lB);
end

%% ======================================================================
function [Gx,Gy] = cell_gradient(W,S)
% Least-squares gradient of every column of W.
Gx = zeros(size(W));  Gy = zeros(size(W));
for k = 1:4
    dW = W(S.nb0(:,k),:) - W;
    Gx = Gx + S.cx(:,k).*dW;
    Gy = Gy + S.cy(:,k).*dW;
end
end

%% ======================================================================
function phi = limiter(r,type)
switch type
    case 'minmod',    phi = max(0,min(1,r));
    case 'vanleer',   phi = (r + abs(r))./(1 + abs(r));
    case 'mc',        phi = max(0,min(min(2*r,0.5*(1+r)),2));
    case 'vanalbada', phi = max(0,(r.^2 + r)./(r.^2 + 1));
    case 'superbee',  phi = max(0,max(min(2*r,1),min(r,2)));
    otherwise, error('Unknown limiter: %s',type);
end
end

%% ======================================================================
function F = physflux(r,u,v,p,nx,ny,gam)
un = u.*nx + v.*ny;
E  = p/(gam-1) + 0.5*r.*(u.^2 + v.^2);
F  = [r.*un, r.*u.*un + p.*nx, r.*v.*un + p.*ny, un.*(E + p)];
end

%% ======================================================================
function F = numflux(WL,WR,nx,ny,gam,type)
% Riemann flux through faces with unit normal (nx,ny), computed in the
% normal/tangential frame and rotated back.
rL = WL(:,1);  pL = WL(:,4);   rR = WR(:,1);  pR = WR(:,4);
unL = WL(:,2).*nx + WL(:,3).*ny;   utL = -WL(:,2).*ny + WL(:,3).*nx;
unR = WR(:,2).*nx + WR(:,3).*ny;   utR = -WR(:,2).*ny + WR(:,3).*nx;

EL = pL/(gam-1) + 0.5*rL.*(unL.^2 + utL.^2);
ER = pR/(gam-1) + 0.5*rR.*(unR.^2 + utR.^2);
UL = [rL, rL.*unL, rL.*utL, EL];
UR = [rR, rR.*unR, rR.*utR, ER];
FL = [rL.*unL, rL.*unL.^2 + pL, rL.*unL.*utL, unL.*(EL + pL)];
FR = [rR.*unR, rR.*unR.^2 + pR, rR.*unR.*utR, unR.*(ER + pR)];
aL = sqrt(gam*pL./rL);   aR = sqrt(gam*pR./rR);

if strcmp(type,'rusanov')
    smax = max(abs(unL) + aL, abs(unR) + aR);
    Fn = 0.5*(FL + FR) - 0.5*smax.*(UR - UL);
else
    % Roe averages -> Einfeldt wave speed estimates
    sq1 = sqrt(rL);  sq2 = sqrt(rR);
    w1 = sq1./(sq1 + sq2);  w2 = 1 - w1;
    HT = w1.*(EL + pL)./rL + w2.*(ER + pR)./rR;
    unT = w1.*unL + w2.*unR;   utT = w1.*utL + w2.*utR;
    aT = sqrt(max((gam-1)*(HT - 0.5*(unT.^2 + utT.^2)),1e-12));
    SL = min(unL - aL, unT - aT);
    SR = max(unR + aR, unT + aT);

    if strcmp(type,'hll')
        Fh = (SR.*FL - SL.*FR + (SL.*SR).*(UR - UL))./(SR - SL);
        m1 = double(SL >= 0);  m3 = double(SR <= 0);  m2 = 1 - m1 - m3;
        Fn = m1.*FL + m2.*Fh + m3.*FR;
    else   % hllc
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
end
F = [Fn(:,1), Fn(:,2).*nx - Fn(:,3).*ny, Fn(:,2).*ny + Fn(:,3).*nx, Fn(:,4)];
end

%% ======================================================================
function [fld,ttl] = field_of(name,W,Gx,Gy,P,Lref,U_inf)
rho = W(:,1);  u = W(:,2);  v = W(:,3);  p = W(:,4);
Vm = hypot(u,v);
switch lower(name)
    case 'mach',      fld = Vm./sqrt(P.gam*p./rho);   ttl = 'Mach number';
    case 'pressure',  fld = p;                        ttl = 'Pressure [Pa]';
    case 'density',   fld = rho;                      ttl = 'Density [kg/m^3]';
    case 'velocity',  fld = Vm;                       ttl = 'Velocity magnitude [m/s]';
    case 'vorticity', fld = (Gx(:,3) - Gy(:,2))*Lref/U_inf;  ttl = 'Vorticity \omega L_{ref}/U_\infty';
    case 'schlieren'
        gm = hypot(Gx(:,1),Gy(:,1));
        fld = exp(-15*gm/max(gm));                    ttl = 'Numerical schlieren';
    case 'entropy',   fld = p./P.p_inf.*(P.rho_inf./rho).^P.gam - 1;  ttl = 'Entropy error (s/s_\infty - 1)';
    otherwise, error('Unknown plotField: %s',name);
end
end
