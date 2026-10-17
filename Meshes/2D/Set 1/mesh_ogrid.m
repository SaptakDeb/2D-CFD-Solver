function m = mesh_ogrid(cfg, body)
% MESH_OGRID  Body-fitted O-grid around a convex outline (blunt bodies).
%   body.P      fine closed polyline, counter-clockwise, starting at the rear point
%   body.brk    indices of forced nodes (start, corners, end)
%   body.hot    clustering hot spots [x y amp width/Lref]
%   body.xb/yb  body centre,  body.Lref reference length,  body.name

mc = cfg.mesh;  dm = cfg.domain;  L = body.Lref;
ds0 = mc.dsWall*L;

mon = struct('aCurv',mc.curvWeight,'kmax',100,'Lref',L,'hot',body.hot);
[xs,ys] = surface_nodes(body.P, body.brk, mc.nSurface, mon);
xs(end) = [];  ys(end) = [];                 % closed curve: drop duplicate

% Outward wall normals (central differences, counter-clockwise ordering)
tx = circshift(xs,-1) - circshift(xs,1);
ty = circshift(ys,-1) - circshift(ys,1);
tl = hypot(tx,ty);  tx = tx./tl;  ty = ty./tl;
nx = ty;   ny = -tx;

% Outer boundary: ellipse covering the requested domain. Wall node i is
% joined to the outer point with the same index fraction around the loop.
xmin = body.xb - dm.up*L;    xmax = body.xb + dm.down*L;
xo = 0.5*(xmin + xmax);      ao = 0.5*(xmax - xmin);
bo = dm.lateral*L;
% Outer angle is spread uniformly with the node index (not with the polar
% angle), so nodes clustered at corners/rear fan out instead of piling up.
th = 2*pi*(0:numel(xs)-1)'/numel(xs);
Pox = xo + ao*cos(th);
Poy = body.yb + bo*sin(th);

[X,Y] = mesh_rays(xs,ys,nx,ny,Pox,Poy,ds0,mc.maxGrowth,mc.blendDist*L);

isWall = false(size(X));  isWall(:,1) = true;
id = reshape(1:numel(X), size(X));
m  = mesh_assemble(X,Y,id,isWall,true);

m.bodyX = xs;  m.bodyY = ys;
m.Lref = L;    m.xb = body.xb;  m.yb = body.yb;   m.name = body.name;
m.topology = 'O-grid';
end
