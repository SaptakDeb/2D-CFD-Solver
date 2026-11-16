%% Test case: NACA 0012, laminar, M = 0.5, Re = 5000, alpha = 0 deg
%  Run from this folder.  Reference: NASA/TM-2016-219003 (see reference_data.txt).
clear; clc; close all;
addpath('..');                         % solver files

P = struct();
P.meshSource = 'file';
P.meshFile   = 'naca0012.msh';         % made by ../gen_naca0012_gmsh.py (chord 1, LE at 0,0)
P.Mach       = 0.5;
P.Re         = 5000;
P.aoa_deg    = 0;
P.flux       = 'hllc';
P.order      = 2;
P.limiter    = 'venkat';
P.cflAdapt   = 'smart';
P.maxIter    = 14000;
P.coefWindow = 2000;
P.coefTol    = 2e-3;
P.avgWindow  = 4000;                   % mean / amplitude of Cl, Cd, Cm over the last 4000 iterations
P.printEvery = 1000;
P.plots      = true;

ns_solver(P);
