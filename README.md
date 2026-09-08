# 2D-CFD-Solver
Exploring numerical methods with Matlab in computational fluid dynamics
---
## Overview
This repository documents my journey of learning how computational fluid dynamics problems are solved numerically.
Rather than simply implementing a solver, the goal is to understand the mathematical reasoning, computational techniques, and engineering decisions behind flow simulation.
The project begins with a two-dimensional potential flow solver around a rectangular solid, solved using Gauss-Seidel iteration in MATLAB/GNU Octave, and will gradually evolve toward more advanced numerical methods and computational fluid dynamics.

---

## Project Goals
* Learn the mathematics behind potential flow theory.
* Understand how the Laplace/potential equation is solved numerically.
* Explore iterative solvers (Gauss-Seidel, SOR).
* Investigate numerical stability and convergence behavior.
* Study solid-boundary treatment (Neumann conditions) and far-field boundaries.
* Extend incompressible potential flow toward compressible flow relations.
* Build intuition for CFD solvers before moving to full Navier-Stokes methods.

# Development Log

# Version 0 — Basic Potential Flow Solver Around a Rectangular Solid

## Code
`v0_potential_flow_rectangle.m`

## Features
- Implemented a 2D potential flow solver using Gauss-Seidel iteration with successive relaxation (SOR).
- 160 × 100 computational grid.
- Rectangular solid embedded in the domain with zero-normal-velocity (Neumann) boundary treatment.
- Far-field Dirichlet boundaries set to freestream potential (phi = U_inf * x).
- Velocity field recovered from the potential via finite-difference gradients.
- Compressible post-processing: temperature, density, pressure, and local Mach number computed from the velocity field.
- Live convergence monitoring (iteration, error, elapsed time) printed during the solve.
- Visualization: velocity magnitude, pressure, Mach number, and velocity vector field, all with the solid outline overlaid.

---

## Notes
- Tested with a single rectangular solid centered in the domain.
- Flow accelerated around the solid as expected, with velocity increasing near the corners.
- Solver converged well below tolerance (1e-7) within the iteration budget at omega = 0.8.
- Compressible property calculations are a post-processing approximation layered on top of an incompressible potential solution, so results are only physically meaningful at low Mach number.
