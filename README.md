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

# Version 2: 2D Euler Solver Around a Rotatable Rectangular Solid (Rusanov FVM)

## Features
- Solves the 2D Euler equations in conservative form, dU/dt + dF/dx + dG/dy = 0 with U = [rho, rho·u, rho·v, E], so density, momentum and energy are all evolved directly instead of through a potential.
- Finite-volume discretisation with a Rusanov (local Lax-Friedrichs) flux, using the maximum local wave speed |u_n| + a at each face.
- Optional second-order MUSCL reconstruction of the primitive variables with a minmod limiter (order = 2). The first 300 iterations run first order for a robust start-up.
- SSP-RK2 time integration with local time stepping, so the solution marches to a steady state faster. The CFL ramps up over the first 300 iterations.
- Much larger domain than Version 1: x from -13.7 to 26.3 and y from -15 to 15, with a longer downstream length for the wake.
- Non-uniform Cartesian grid: a uniform refined box (h0 = 0.025) around the body, then geometric stretching (ratio 1.06) out to the boundaries, about 211 × 202 cells at default settings.
- Adjustable angle of attack (aoa_deg). The rectangle is rotated about its own centre, positive nose-up, while the freestream stays along +x.
- Rotated geometry is represented as a staircase of solid cells. Solid faces act as slip walls by mirroring the normal velocity of the neighbouring fluid cell.
- Far-field boundaries: freestream at the inlet, top and bottom; zero gradient at the outlet with the pressure fixed at p_inf.
- Density and pressure floors (1% of freestream) guard against negative values, and the reconstruction falls back to first order wherever it would give negative density or pressure.
- Lift and drag coefficients from the pressure force on the body surface, based on the body length Lb, and printed live with the residual.
- Relative density residual (divided by its value at the first iteration) used as the convergence measure.
- Live progress (iteration, residual, Cl, Cd, elapsed time) printed during the solve.
- Visualization: velocity magnitude, pressure, Mach number, velocity vector field, and a convergence and force history plot, zoomed on the body with the solid outline overlaid.

---

## Notes
- Tested in GNU Octave 8.4 on a coarser grid (h0 = 0.05) for short runs, and the solver stayed stable. The full default grid and the plotting section were checked less thoroughly.
- The solid is a staircase on a Cartesian grid, so Cl and Cd are approximate, and the extra corners it creates affect the wake.
- The sharp corners of the rectangle act as vortex-sheet sources, so the wake and recirculation bubble depend on the numerical dissipation of the scheme and change with h0 and order. Treat the wake size and shape as scheme-dependent, not as a grid-converged physical result.
- The residual usually plateaus (around 1e-3 in the coarse tests) instead of reaching the tolerance. The run then stops at maxIter, and Cl and Cd may keep oscillating. In that case the "steady" result is one snapshot of an unsteady wake.
- At AoA = 0 the grid, geometry and scheme are symmetric, so the solution stays symmetric: Cl is close to zero and the wake forms a twin-vortex bubble behind the body. That branch may be unstable physically. A small asymmetry (for example aoa_deg = 0.1) is what would trigger alternating shedding.
- The model is still inviscid, so there is no boundary layer, no viscous separation and no viscous drag. The nonzero Cd on this bluff body comes from the vortex-sheet wake and numerical dissipation, and it is not skin friction.
- A smooth body (ellipse or airfoil) at low Mach number should give attached flow with Cd close to zero, so it is a good validation case. This is planned for later.
- The speed of the run is set by the fine cells. Raise h0 or set order = 1 for quick previews.
- The plotted vector field samples the stretched grid, so the columns downstream look unevenly spaced. That is a plotting effect and not a flow feature.
- If the solver diverges or oscillates, lower CFL first (0.5 was the value used in testing; the current default is 0.8), then increase nFirst or switch to order = 1.
---
