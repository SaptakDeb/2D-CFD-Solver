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



# Version 1: Full Potential Solver Around a Rectangular Solid (Red-Black SOR)

## Features
- Solves the 2D full potential equation in conservative form, d(rho·phi_x)/dx + d(rho·phi_y)/dy = 0, so compressibility is part of the solve rather than added afterward.
- Finite-volume discretisation on the 160 × 100 grid, with face densities averaged from neighbouring nodes.
- Red-black SOR iteration (omega = 1.85), vectorised so it runs much faster in Octave than nested loops.
- Isentropic density update each iteration, under-relaxed (omega_rho = 0.8). The speed is capped at about Mach 0.85 inside the density law, and the number of capped nodes is reported.
- Rectangular solid with a zero-flux face wherever a fluid node neighbours the solid.
- Inlet and outlet fixed at the freestream potential (phi = U_inf · x). Top and bottom boundaries use zero normal velocity (Neumann).
- Velocity recovered with the same ghost-value treatment next to the solid, so values inside the body no longer contaminate the surface.
- Scaled convergence error (divided by U_inf · dx), so the tolerance does not depend on grid size or speed.
- Conservation check comparing mass flux through the inlet and the outlet.
- Live convergence monitoring (iteration, error, elapsed time) printed during the solve.
- Visualization: velocity magnitude, pressure, Mach number, velocity vector field, and convergence history, with the solid outline overlaid.

---

## Notes
- Tested with a single rectangular solid in the domain, run in GNU Octave.
- Solver converged to the scaled tolerance (1e-6) within the iteration budget at omega = 1.85 and omega_rho = 0.8.
- Flow accelerated around the solid as expected, with the highest velocities at the corners.
- The sharp corners of the rectangle produce very large discrete velocities. Treat the corner values as a numerical artifact, not as physical results.
- The inlet and outlet mass fluxes agreed closely, which supports the conservative formulation.
- The model is still inviscid and irrotational, so there is no wake or separation. It also cannot capture shocks, which is why the Euler equations are the planned Version 2.
- If the solver diverges or oscillates for other geometries or Mach numbers, lower omega to about 1.7 first, then lower omega_rho.
---

