function m = mesh_cylinder(cfg)
% MESH_CYLINDER  O-grid around a circular cylinder (Lref = diameter).
R  = cfg.body.radius;   xb = cfg.body.xc;   yb = cfg.body.yc;

th = linspace(0,2*pi,6000)';
P  = [xb + R*cos(th), yb + R*sin(th)];
P(end,:) = P(1,:);

body.P    = P;
body.brk  = [1 size(P,1)];
body.hot  = [xb+R, yb, cfg.mesh.wakeBoost, 0.3];      % rear point -> wake side
body.xb = xb;  body.yb = yb;  body.Lref = 2*R;  body.name = 'cylinder';
m = mesh_ogrid(cfg, body);
end
