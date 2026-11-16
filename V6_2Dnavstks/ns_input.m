clear; clc; close all;

%% ========================================================================
%  INPUT FILE  -  2D compressible Navier-Stokes / Euler, unstructured FVM
%
%  Files (keep them in one folder, then run this file):
%    ns_input.m        this file: all settings
%    ns_solver.m       solver (flux, reconstruction, time stepping, CFL adaptation, forces)
%    read_mesh.m       mesh reader: Gmsh .msh (2.2 / 4.x ASCII), mesh2d.py text files,
%                      or built-in O-grids
%    build_fvmesh.m    edges, ghost cells, boundary classification
%    mesh_*.m, ogrid_extrude.m   built-in O-grid generators (cylinder, ellipse, naca0012, rectangle)
%    gen_naca0012_gmsh.py        example: NACA 0012 mesh with boundary layer (needs gmsh)
%  ========================================================================

P = struct();

%% Mesh
P.meshSource = 'file';                 % 'file' or 'builtin'
P.meshFile   = 'testcase_naca0012/naca0012.msh';   % .msh (Gmsh ASCII) or the base name of the mesh2d.py text files
P.wallTags   = [];                     % physical tags of the body (empty = detect from group names)
P.farTags    = [];                     % physical tags of the far field
P.Lref       = [];                     % reference length (empty = chord of the wall nodes)
P.momentRef  = [];                     % moment reference point [x y] (empty = quarter chord)

% only for P.meshSource = 'builtin'
P.geometry   = 'cylinder';             % 'cylinder', 'ellipse', 'naca0012', 'rectangle'
P.Rc = 0.5;  P.Lb = 1.0;  P.Hb = 0.7;  P.naca_t = 0.12;
P.Ni = 160;  P.Nj = 90;  P.h1_Lref = 0.015;  P.Rfar_Lref = 20;  P.blendPower = 1.5;

%% Flow
P.Mach       = 0.5;
P.aoa_deg    = 0;                      % rotates the freestream, the mesh stays fixed
P.Re         = 5000;                   % Inf = Euler equations (slip wall)
P.Pr         = 0.72;
P.gamma      = 1.4;
P.viscosity  = 'constant';             % 'constant' or 'sutherland'
P.Tinf_dim   = 288.15;                 % [K], only used by Sutherland
P.wall       = 'adiabatic';            % 'adiabatic' or 'isothermal' (viscous runs)
P.Tw_ratio   = 1.0;                    % Tw / T_inf for the isothermal wall

%% Numerical methods
P.flux       = 'hllc';                 % 'rusanov' or 'hllc'
P.order      = 2;                      % 1 = first order, 2 = linear reconstruction
P.limiter    = 'venkat';               % 'venkat' or 'barth'
P.venK       = 3;                      % Venkatakrishnan constant (larger = less limiting)
P.timeScheme = 'rk3';                  % 'rk2' or 'rk3'
P.nFirst     = 200;                    % first-order start-up iterations

%% CFL control
P.cflAdapt     = 'smart';              % 'off'   fixed CFL with a start-up ramp
                                       % 'smart' grow the CFL while steps are accepted, cut it and
                                       %         roll back when a step gives negative density /
                                       %         pressure, NaN, or a residual jump
                                       % 'ser'   as 'smart', CFL also follows the residual drop
P.CFL          = 1.0;                  % start value (the fixed value for 'off')
P.cflMin       = 0.1;
P.cflMax       = 3.0;
P.cflUp        = 1.05;                 % growth factor
P.cflDown      = 0.5;                  % reduction factor on a rejected step
P.cflGrowEvery = 10;                   % accepted steps between increases ('smart')
P.cflJump      = 5;                    % reject a step if the residual grows by more than this factor
P.cflGuardIter = 50;                   % residual-jump test starts after this many iterations
P.cflSerExp    = 1.0;                  % exponent of the residual ratio ('ser')
P.rampIter     = 300;                  % start-up ramp length ('off')

%% Run
P.mode       = 'steady';               % 'steady' (local time step) or 'transient' (global time step)
P.maxIter    = 30000;
P.tol        = 1e-8;                   % steady: density residual relative to its largest value
P.coefWindow = 1000;                   % steady: also stop when Cl, Cd, Cm vary by less than coefTol
P.coefTol    = 2e-3;                   %         over the last coefWindow iterations
P.printEvery = 200;
P.plots      = true;

% transient
P.tEnd_conv  = 100;                    % end time, t* = t U_inf / Lref
P.plotEvery  = 500;
P.plotField  = 'vorticity';            % 'vorticity','mach','pressure','density','velocity','temperature'
P.saveFrames = false;
P.perturbAmp = 1e-2;                   % start-up wake perturbation (fraction of U_inf)

%% Run the solver (results are also left in the workspace variable ns_results)
ns_solver(P);
