clear; clc; close all;
here = fileparts(mfilename('fullpath'));
if ~isempty(here), addpath(here); end

%% Flow conditions
cfg.flow.gamma   = 1.4;
cfg.flow.R       = 287;
cfg.flow.M_inf   = 0.3;
cfg.flow.T_inf   = 300;
cfg.flow.p_inf   = 101325;
cfg.flow.aoa_deg = 0;               % angle of attack [deg], positive = nose up
                                    % (the flow is rotated, the mesh stays fixed)

%% Body
cfg.shape = 'naca0012';             % 'naca0012', 'cylinder', 'ellipse', 'rectangle'

cfg.body.xc        = 1.3;           % body centre x
cfg.body.yc        = 0.0;           % body centre y
cfg.body.length    = 1.0;           % chord (naca0012) / length (ellipse, rectangle)  = Lref
cfg.body.height    = 0.7;           % thickness (ellipse, rectangle)
cfg.body.radius    = 0.5;           % cylinder (Lref = diameter)
cfg.body.tc        = 0.12;          % naca0012 thickness ratio

%% Domain (in units of Lref, measured from the body centre)
cfg.domain.up      = 15;            % upstream
cfg.domain.down    = 30;            % downstream (long, for the wake)
cfg.domain.lateral = 15;            % half height

%% Mesh
cfg.mesh.nSurface   = 240;          % intervals around the body
cfg.mesh.dsWall     = 0.02;         % first cell height / Lref (inviscid: no boundary layer)
cfg.mesh.maxGrowth  = 1.08;         % wall-normal growth ratio (sets the number of layers)
cfg.mesh.blendDist  = 2;            % grid lines follow the wall normal over this distance / Lref
cfg.mesh.wakeGrowth = 1.05;         % streamwise growth in the wake (C-grid)
cfg.mesh.curvWeight = 0.7;          % clustering with curvature (nose, tail)
cfg.mesh.edgeBoost  = 3;            % extra clustering at trailing edge / corners
cfg.mesh.wakeBoost  = 2;            % extra clustering on the wake side (O-grids)
cfg.mesh.plot       = true;         % plot the mesh before solving
cfg.mesh.meshOnly   = false;        % true = stop after the mesh

%% Numerical methods
cfg.num.mode       = 'steady';      % 'steady' (local time step) or 'transient' (global time step)
cfg.num.flux       = 'hllc';        % 'rusanov', 'hll', 'hllc'
cfg.num.order      = 2;             % 1 = first order, 2 = MUSCL
cfg.num.limiter    = 'vanleer';     % 'minmod', 'vanleer', 'mc', 'vanalbada', 'superbee'
cfg.num.timeScheme = 'rk3';         % 'rk2', 'rk3'
cfg.num.farfield   = 'characteristic';   % 'characteristic' or 'freestream'

cfg.num.CFL        = 0.8;
cfg.num.nFirst     = 200;           % first-order start-up iterations
cfg.num.rampIter   = 300;           % CFL ramp length
cfg.num.maxIter    = 6000;
cfg.num.tol        = 1e-5;          % steady only: relative density residual
cfg.num.printEvery = 10;

%% Transient settings
cfg.num.tEnd_conv  = 30;            % end time in convective units  t* = t U_inf / Lref
cfg.num.plotEvery  = 20;            % plot every N iterations
cfg.num.plotField  = 'vorticity';   % 'mach','pressure','density','velocity','vorticity','schlieren','entropy'
cfg.num.saveFrames = false;         % true = save frame_XXXXX.png
cfg.num.perturbAmp = 1e-2;          % start-up wake perturbation (fraction of U_inf)

%% Mesh
mesh = feval(['mesh_' lower(cfg.shape)], cfg);
fprintf('Mesh (%s, %s): %d cells, cell area %.2e .. %.2e\n', ...
        mesh.name, mesh.topology, numel(mesh.vol), min(mesh.vol), max(mesh.vol));

if cfg.mesh.plot
    figure('Name','Mesh');
    patch('Faces',mesh.cells,'Vertices',[mesh.xn mesh.yn],'FaceColor','none','EdgeColor',[0.3 0.3 0.3]);
    hold on;  line(mesh.wx,mesh.wy,'Color','r','LineWidth',1.5);
    axis equal;  xlabel('x');  ylabel('y');  title('Mesh (zoom with the axes tools)');
    xlim([mesh.xb - 2*mesh.Lref, mesh.xb + 4*mesh.Lref]);
    ylim([mesh.yb - 2*mesh.Lref, mesh.yb + 2*mesh.Lref]);
    drawnow;
end
if cfg.mesh.meshOnly, return; end

%% Solve and post-process
sol = euler_solver(cfg, mesh);
euler_post(cfg, mesh, sol);
