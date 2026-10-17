function m = mesh_rectangle(cfg)
% MESH_RECTANGLE  O-grid around a sharp-cornered rectangle (Lref = length).
%   The four corners are forced nodes and carry clustering hot spots, so
%   they are resolved exactly instead of being cut by a staircase mask.
a  = 0.5*cfg.body.length;   b = 0.5*cfg.body.height;
xb = cfg.body.xc;           yb = cfg.body.yc;

c = [xb+a yb; xb+a yb+b; xb-a yb+b; xb-a yb-b; xb+a yb-b; xb+a yb];   % CCW, from rear mid-point
n = 1500;
P = [xb+a, yb];  brk = 1;
for k = 1:5
    seg = [linspace(c(k,1),c(k+1,1),n)', linspace(c(k,2),c(k+1,2),n)'];
    P   = [P; seg(2:end,:)];
    brk(end+1) = size(P,1);
end

bm = cfg.mesh.edgeBoost;
body.P    = P;
body.brk  = brk;
body.hot  = [xb+a yb+b bm 0.05; xb-a yb+b bm 0.05; xb-a yb-b bm 0.05; xb+a yb-b bm 0.05; ...
             xb+a yb cfg.mesh.wakeBoost 0.3];
body.xb = xb;  body.yb = yb;  body.Lref = cfg.body.length;  body.name = 'rectangle';
m = mesh_ogrid(cfg, body);
end
