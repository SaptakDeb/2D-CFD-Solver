function [X,Y] = mesh_rays(xs,ys,nx,ny,Pox,Poy,ds0,growth,delta)
% MESH_RAYS  Body-fitted grid lines from the wall to the outer boundary.
%   Each wall node (xs,ys) is joined to its outer boundary node (Pox,Poy).
%   Near the wall the line leaves along the wall normal (nx,ny) so cells
%   are orthogonal; over the distance delta it blends smoothly into the
%   straight line to the outer node. Node spacing along every line is
%   geometric with first cell height ds0.
%   X, Y are (number of wall nodes) x (NJ+1).

ex = Pox - xs;  ey = Poy - ys;
Lr = hypot(ex,ey);
ex = ex./Lr;    ey = ey./Lr;

NJ = ceil(log(1 + max(Lr)*(growth-1)/ds0)/log(growth));
s  = geo_ratio(Lr, NJ, ds0);
d  = Lr.*s;                                 % distance from the wall

t  = min(d/delta, 1);
w  = t.^2.*(3 - 2*t);                       % smoothstep, 0 at wall, 1 beyond delta
vx = (1-w).*nx + w.*ex;
vy = (1-w).*ny + w.*ey;
vl = hypot(vx,vy);

X = xs + d.*vx./vl;
Y = ys + d.*vy./vl;
end
