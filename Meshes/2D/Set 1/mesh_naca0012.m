function m = mesh_naca0012(cfg)
% MESH_NACA0012  C-grid around a NACA 0012 with a closed (sharp) trailing edge.
%   Layout of the i-index: lower wake (far -> TE), lower surface, upper
%   surface, upper wake (TE -> far). The wake cut is closed by giving the
%   two wake rows the same node ids, so the solver sees an ordinary
%   interior boundary there.

mc = cfg.mesh;  dm = cfg.domain;
L  = cfg.body.length;   xb = cfg.body.xc;   yb = cfg.body.yc;
tt = cfg.body.tc;
nS  = mc.nSurface;                 % intervals along the whole surface
ds0 = mc.dsWall*L;

% Outline: lower TE -> LE -> upper TE (closed-TE coefficient -0.1036)
xn = 0.5*(1 - cos(linspace(0,pi,4000)))';
yt = 5*tt*(0.2969*sqrt(xn) - 0.1260*xn - 0.3516*xn.^2 + 0.2843*xn.^3 - 0.1036*xn.^4);
P  = [flipud([xn -yt]); [xn(2:end) yt(2:end)]];
P  = [xb + L*(P(:,1) - 0.5), yb + L*P(:,2)];
xTE = xb + 0.5*L;

mon = struct('aCurv',mc.curvWeight,'kmax',100,'Lref',L, ...
             'hot',[xTE yb mc.edgeBoost 0.04]);
[xs,ys] = surface_nodes(P,[1 size(P,1)],nS,mon);
ns = numel(xs);

% Outward normals (this ordering is clockwise: n = (-ty, tx))
tx = zeros(ns,1);  ty = zeros(ns,1);
tx(2:end-1) = xs(3:end) - xs(1:end-2);   ty(2:end-1) = ys(3:end) - ys(1:end-2);
tx(1) = xs(2) - xs(1);   ty(1) = ys(2) - ys(1);
tx(end) = xs(end) - xs(end-1);   ty(end) = ys(end) - ys(end-1);
tl = hypot(tx,ty);  tx = tx./tl;  ty = ty./tl;
snx = -ty;   sny = tx;

% Near the trailing edge the wall normals lean towards the wake (by the
% wedge half angle) and would cross the vertical wake columns. Blend them
% to vertical there so the lines fan out cleanly from the TE.
[~,iLE] = min(xs);
vs = [-ones(iLE,1); ones(ns-iLE,1)];
cT = exp(-((xTE - xs)/(0.15*L)).^2);
snx = (1-cT).*snx;
sny = (1-cT).*sny + cT.*vs;
sl  = hypot(snx,sny);   snx = snx./sl;   sny = sny./sl;

% Wake nodes: geometric growth starting from the trailing-edge spacing
xmax = xb + dm.down*L;
h0 = hypot(xs(2) - xs(1), ys(2) - ys(1));
r  = mc.wakeGrowth;
LW = xmax - xTE;
nW = ceil(log(1 + LW*(r-1)/h0)/log(r));
ww = h0*r.^(0:nW-1)';   ww = ww*(LW/sum(ww));
xw = xTE + cumsum(ww);

% Wall / wake node row (j = 1)
xi  = [flipud(xw); xs; xw];
yi  = [zeros(nW,1); ys; zeros(nW,1)];
nxi = [zeros(nW,1); snx; zeros(nW,1)];
nyi = [-ones(nW,1); sny; ones(nW,1)];
NIn = numel(xi);

% Outer boundary: straight lines above/below the wake, half ellipse in front
H    = dm.lateral*L;
aUp  = xTE - (xb - dm.up*L);
beta = pi*(0:nS)'/nS;
Pox  = [flipud(xw); xTE - aUp*sin(beta); xw];
Poy  = [-H*ones(nW,1); -H*cos(beta); H*ones(nW,1)];

[X,Y] = mesh_rays(xi,yi,nxi,nyi,Pox,Poy,ds0,mc.maxGrowth,mc.blendDist*L);

% Node numbering: lower and upper wake rows (and the TE) are the same nodes
id = reshape(1:numel(X), size(X));
for i = 1:nW+1
    id(NIn+1-i,1) = id(i,1);
end
isWall = false(size(X));
isWall(nW+1:nW+ns,1) = true;

m = mesh_assemble(X,Y,id,isWall,false);
m.bodyX = xs;  m.bodyY = ys;
m.Lref = L;    m.xb = xb;  m.yb = yb;   m.name = 'naca0012';
m.topology = 'C-grid';
end
