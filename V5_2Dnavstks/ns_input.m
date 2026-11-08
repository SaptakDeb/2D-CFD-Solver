clear; clc; close all;

%% ========================================================================
%  INPUT FILE  -  2D laminar compressible Navier-Stokes (Octave / MATLAB)
%
%  Files (keep them in the same folder, then run this file):
%    ns_input.m        this file: all user settings
%    ns_solver.m       finite-volume solver
%    ogrid_extrude.m   builds the body-fitted O-grid from a body curve
%    mesh_cylinder.m   mesh file, cylinder
%    mesh_ellipse.m    mesh file, ellipse
%    mesh_naca0012.m   mesh file, NACA 00xx airfoil
%    mesh_rectangle.m  mesh file, rectangle
%  ========================================================================

P = struct();

%% Geometry
P.geometry = 'cylinder';     % 'cylinder', 'ellipse', 'naca0012', 'rectangle'
P.Rc       = 0.5;            % cylinder radius
P.Lb       = 1.0;            % ellipse / rectangle length, NACA chord
P.Hb       = 0.7;            % ellipse / rectangle thickness
P.naca_t   = 0.12;           % NACA thickness ratio
P.aoa_deg  = 0;              % angle of attack [deg], positive = nose up (rotates the geometry)

%% Flow
P.Mach       = 0.3;
P.Re         = 40;           % Re = rho U Lref / mu  (Lref = diameter, length or chord)
P.Pr         = 0.72;
P.gamma      = 1.4;
P.viscosity  = 'constant';   % 'constant' or 'sutherland'
P.Tinf_dim   = 288.15;       % freestream temperature [K], only used by Sutherland
P.wall       = 'adiabatic';  % 'adiabatic' or 'isothermal'
P.Tw_ratio   = 1.0;          % wall temperature / freestream temperature (isothermal only)

%% Mesh (O-grid: i around the body, j from the wall to the far field)
P.Ni         = 160;          % cells around the body (even number, needed by the NACA mesh)
P.Nj         = 90;           % cells from wall to far field
P.h1_Lref    = 0.015;        % first cell height at the wall / Lref (use ~Re^-0.5 / 10 or smaller)
P.Rfar_Lref  = 20;           % far-field radius / Lref
P.blendPower = 1.5;          % grid lines follow the wall normal longer for larger values

%% Numerical methods
P.flux       = 'hllc';       % 'rusanov' or 'hllc'
P.order      = 2;            % 1 = first order, 2 = MUSCL
P.limiter    = 'vanleer';    % 'minmod', 'vanleer', 'mc', 'vanalbada'
P.timeScheme = 'rk3';        % 'rk2' or 'rk3'
P.CFL        = 0.8;
P.nFirst     = 200;          % first-order start-up iterations
P.rampIter   = 300;          % CFL ramp length

%% Run
P.mode       = 'steady';     % 'steady' (local time step) or 'transient' (global time step)
P.maxIter    = 30000;
P.tol        = 1e-6;         % steady: relative density residual
P.printEvery = 200;

% Transient / vortex shedding (Re ~ 100 for the cylinder, mode = 'transient')
P.tEnd_conv  = 100;          % end time in t* = t U_inf / Lref
P.plotEvery  = 500;          % plot every N iterations
P.plotField  = 'vorticity';  % 'vorticity','mach','pressure','density','velocity','temperature'
P.saveFrames = false;        % true = save frame_XXXXX.png
P.perturbAmp = 1e-2;         % start-up wake perturbation (fraction of U_inf), 0 = none

%% Run the solver
ns_solver(P);
