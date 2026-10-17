function m = mesh_ellipse(cfg)
% MESH_ELLIPSE  O-grid around an ellipse (length x height, Lref = length).
a  = 0.5*cfg.body.length;   b = 0.5*cfg.body.height;
xb = cfg.body.xc;           yb = cfg.body.yc;

th = linspace(0,2*pi,6000)';
P  = [xb + a*cos(th), yb + b*sin(th)];
P(end,:) = P(1,:);

body.P    = P;
body.brk  = [1 size(P,1)];
body.hot  = [xb+a, yb, cfg.mesh.wakeBoost, 0.3];
body.xb = xb;  body.yb = yb;  body.Lref = cfg.body.length;  body.name = 'ellipse';
m = mesh_ogrid(cfg, body);
end
