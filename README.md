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
# Version 3: 2D Compressible Euler Solver (HLLC FVM)

## Features
- Solves the 2D compressible Euler equations using a finite-volume method.
- Selectable Rusanov, HLL and HLLC flux solvers.
- First- and second-order MUSCL reconstruction with selectable slope limiters.
- SSP-RK2/RK3 time integration with configurable CFL.
- Supports rectangle, cylinder, ellipse and NACA 0012 geometries.
- Adjustable angle of attack and body dimensions.
- Non-uniform Cartesian grid with refined body and wake regions.
- Approximate slip-wall and characteristic far-field boundary conditions.
- Calculates lift, drag, Mach number, pressure coefficient, vorticity and entropy error.
- Transient flow visualisation with optional frame saving.
- Estimates vortex-shedding frequency and Strouhal number using Fourier analysis.
- Generates flow-field plots and convergence histories.

---

## Notes
- Inviscid Euler model: no physical viscosity, boundary layers, skin friction or viscous separation.
- Cartesian staircase geometry makes surface forces and near-wall flow grid-dependent.
- Wake development and vortex shedding depend on numerical dissipation, resolution and initial perturbations.
- Strouhal-number estimates depend on simulation duration and spectral resolution.
- Grid convergence and validation are required before treating results as quantitatively reliable.
---
# Version 4: 2D Compressible Euler Solver (HLLC FVM, body-fitted mesh)

Version 4 splits the code into separate file groups: an input file, a solver, and separate mesh generators (one per geometry). The solver no longer builds its own grid. It reads a mesh produced by a `mesh_*.m` file, so geometries and meshes can change without touching the solver.


## Features

* Solves the 2D compressible Euler equations using a cell-centred finite-volume method on a body-fitted quad mesh.
* Selectable Rusanov, HLL and HLLC flux solvers.
* First- and second-order MUSCL reconstruction from least-squares gradients, with selectable limiters (minmod, van Leer, MC, van Albada, superbee).
* SSP-RK2/RK3 time integration with configurable CFL; steady (local time step) or transient (global time step).
* Supports rectangle, cylinder, ellipse and NACA 0012 geometries.
* Adjustable angle of attack (the flow is rotated; the mesh stays fixed) and body dimensions.
* Body-fitted meshes: C-grid for the NACA 0012, O-grids for the other shapes.
* Curvature- and hot-spot-based node clustering (nose, trailing edge, corners, wake side); wall-normal layers set by first cell height and growth ratio.
* Slip wall applied as a pressure flux on the exact body surface; characteristic (or freestream) far field.
* Calculates lift, drag, Mach number, pressure coefficient, surface Cp, vorticity and entropy error.
* Transient flow visualisation with optional frame saving.
* Estimates vortex-shedding frequency and Strouhal number using Fourier analysis.
* Generates flow-field plots and convergence histories.

## Notes

* Inviscid Euler model: no physical viscosity, boundary layers, skin friction or viscous separation.
* Body-fitted meshes remove the staircase geometry of version 3, so surface forces are taken on the true surface. Results still depend on the mesh.
* Mesh refinement is geometry-driven, not solution-adaptive: it does not refine from the flow field.
* O-grid wake resolution (cylinder, ellipse, rectangle) is coarser than in the C-grid NACA wake.
* Wake development and vortex shedding depend on numerical dissipation, resolution and initial perturbations.
* Transient runs use explicit time stepping, so the step count is set by the smallest cell. The solver prints an estimate at the start.
* Strouhal-number estimates depend on simulation duration and spectral resolution.
* Grid convergence and validation are required before treating results as quantitatively reliable.

* ## Additional Information

Put all files in one folder, then in MATLAB/Octave run:   `run_euler `

Edit settings only in `run_euler.m`. Set `cfg.mesh.meshOnly = true` to preview the mesh without solving.

## File structure

| Group | File | Role |
|---|---|---|
| Input | `run_euler.m` | All user settings; builds the mesh, runs the solver, plots |
| Solver | `euler_solver.m` | Finite-volume Euler solver (fluxes, reconstruction, time marching, boundaries, forces) |
| Post-processing | `euler_post.m` | Final plots |
| Mesh (geometry) | `mesh_naca0012.m` | C-grid, sharp trailing edge, wake cut |
| | `mesh_cylinder.m`, `mesh_ellipse.m`, `mesh_rectangle.m` | O-grids |
| Mesh (shared tools) | `mesh_ogrid.m` | O-grid builder used by cylinder, ellipse, rectangle |
| | `mesh_assemble.m` | Converts the node array into cells, faces, normals, neighbours and gradient coefficients for the solver |

---


# Version 5: 2D Compressible Navier-Stokes Solver (HLLC FVM, body-fitted O-grid)

Version 5 keeps the file-group layout of Version 4 (an input file, a solver, and separate mesh generators, one per geometry) and adds viscous terms. The solver does not build its own grid. It reads a mesh produced by a `mesh_*.m` file, so geometries and meshes can change without touching the solver.

## Features
- Solves the 2D laminar compressible Navier-Stokes equations using a cell-centred finite-volume method on a body-fitted quad O-grid (i around the body, j from the wall to the far field).
- Selectable Rusanov and HLLC convective fluxes.
- First- and second-order MUSCL reconstruction (in index space) with selectable limiters (minmod, van Leer, MC, van Albada). A first-order fallback is used wherever the reconstruction would give negative density or pressure.
- Viscous fluxes from Green-Gauss cell gradients, with face gradients corrected along the cell-to-cell direction to avoid odd-even decoupling.
- Constant or Sutherland viscosity, constant Prandtl number, and Reynolds number set from the body reference length.
- No-slip wall, adiabatic or isothermal, with the wall pressure flux applied on the exact body surface.
- Characteristic (Riemann invariant) far-field boundary condition.
- SSP-RK2/RK3 time integration with configurable CFL; steady (local time step, including a viscous limit) or transient (global time step). Start-up uses a first-order phase and a CFL ramp.
- Supports cylinder, ellipse, NACA 00xx (default 0012) and rectangle geometries, each with its own mesh file.
- Adjustable angle of attack (the geometry and mesh are rotated; the freestream stays along +x) and body dimensions.
- Wall-normal layers set by the first cell height and a geometric growth ratio. Surface points are clustered at the NACA leading and trailing edges and at the rectangle corners.
- Calculates lift, drag (split into pressure and friction), Mach number, pressure coefficient, vorticity, wall Cp and wall Cf.
- Estimates the recirculation length on the rear symmetry line at zero angle of attack.
- Transient flow visualisation (vorticity, Mach, pressure, density, velocity, temperature) with optional frame saving.
- Estimates vortex-shedding frequency and Strouhal number using Fourier analysis of the Cl signal and a wake probe.
- Generates flow-field plots, a mesh view, wall distributions and convergence histories.

---

## Notes
- Laminar model only: no turbulence model, so it is only meaningful for Reynolds numbers where the flow is laminar (the cylinder is validated against values around Re = 40 to 200). It is also 2D, so there is no three-dimensional shedding.
- Tested in GNU Octave 8.4 on coarse grids. A steady cylinder at Re = 40 gave a recirculation length of about 2.1 diameters, close to published values, but Cd was still settling and the residual had not converged after 3000 iterations. The NACA 0012, rectangle and ellipse cases at 5 degrees were run for 300 iterations only, and the transient cylinder case was run to t* = 0.5. No long shedding run has been checked.
- Explicit time stepping with a wall-resolved mesh is stiff: the step is set by the smallest wall cell, so transient runs need on the order of a thousand steps per convective time. Steady mode and coarser grids are much cheaper. The solver does not print a step-count estimate.
- All geometries use an O-grid, so the wake is coarse far from the body, and vortex shedding may be damped numerically downstream. A C-grid or wake refinement would be needed for a clean vortex street.
- Straight algebraic grid lines blended toward the wall normal can become skewed at sharp corners (rectangle) and at the NACA trailing edge. The solver stops with an error if a cell has zero area, but it does not check for skewness.
- Wall resolution matters: the first cell height should be much smaller than the boundary layer thickness (roughly Re^-0.5 times 0.1 of the reference length or less).
- Mesh refinement is geometry-driven, not solution-adaptive: it does not refine from the flow field.
- The far-field boundary uses 1D Riemann invariants, which is less accurate if the wake crosses the outer boundary at a steep angle. A larger far-field radius reduces the effect.
- Strouhal-number estimates depend on simulation duration and spectral resolution.
- Grid convergence and validation are required before treating results as quantitatively reliable.

---

# 2D Cartesian Cut-Cell Mesh Generator (Python, Gmsh)

A standalone 2D mesh generator written in Python. It constructs a Cartesian background grid and clips cells intersecting the body geometry, retaining only the fluid region. Boundary-intersecting cells are subdivided into triangles, allowing the body surface to be represented exactly without requiring a body-fitted coordinate system.

The mesher operates independently of the flow solver and exports Gmsh 2.2 meshes alongside plain-text mesh files compatible with GNU Octave.

## Features
- Generates a 2D Cartesian quadrilateral background grid with configurable uniform cell size.
- Supports arbitrary 2D geometries supplied as text files containing x-y coordinate pairs. Header lines are ignored.
- Clips cells intersecting the body boundary against the fluid region, preserving the exact body surface.
- Splits clipped boundary cells into triangles while retaining quadrilateral cells in unaffected regions.
- Supports built-in circle, ellipse, rectangle and NACA 4-digit geometries.
- Allows geometry rotation through an adjustable angle of attack, with positive angles corresponding to nose-up rotation.
- Configurable computational domain limits (`xmin`, `xmax`, `ymin`, `ymax`) and background grid spacing (`h`).
- Exports Gmsh 2.2 `.msh` files with separate physical boundary tags for the body and far-field boundaries.
- Exports plain-text mesh files (`_nodes.txt`, `_tri.txt`, `_quad.txt` and `_bnd.txt`) for direct loading into GNU Octave.
- Includes a separate mesh inspection utility for checking cell connectivity, geometry and mesh topology.
- Detects zero-area and negative-area cells, repeated nodes, duplicate cells, coincident nodes and excessively small cell angles.
- Checks for edges shared by more than two cells, gaps, hanging nodes and inconsistencies in total fluid-domain area.
- Returns exit code 1 when mesh errors are detected.
- Includes a sample NACA 2412 geometry file in Selig coordinate format.

---

# Version 6: 2D Compressible Navier-Stokes Solver (unstructured FVM, Gmsh / mesh2d meshes)

Version 6 reads an external mesh instead of building its own grid, and runs on mixed triangle / quad meshes. Settings live in one input file, the solver is separate from the mesh readers, and the built-in O-grid generators of Version 5 are still available as one mesh file per geometry.

## Files
- `ns_input.m`: all user settings. Run this file.
- `ns_solver.m`: solver (fluxes, reconstruction, viscous terms, boundary conditions, time stepping, CFL control, forces, plots).
- `read_mesh.m`: mesh reader (Gmsh ASCII 2.2 / 4.x, `mesh2d.py` text files, or built-in O-grids).
- `build_fvmesh.m`: cells, edges, ghost cells, boundary classification, gradient weights.
- `mesh_cylinder.m`, `mesh_ellipse.m`, `mesh_naca0012.m`, `mesh_rectangle.m`, `ogrid_extrude.m`: built-in O-grid generators.
- `gen_naca0012_gmsh.py`: NACA 0012 mesh with a boundary layer (needs `pip install gmsh numpy`).
- `testcase_naca0012/`: mesh, two input scripts, reference data, and the results I obtained.

## Features
- Solves the 2D compressible Navier-Stokes equations (laminar) or the Euler equations (`Re = Inf`) with a cell-centred finite-volume method.
- Reads Gmsh ASCII meshes (format 2.2 and 4.x), the text files written by `mesh2d.py`, or the built-in structured O-grids, and converts them all to one edge-based form.
- Triangles and quads in the same mesh.
- Boundary classification from the physical-group names (wall, body, airfoil, farfield, inlet, outlet ...), from `P.wallTags` / `P.farTags`, or from position if the mesh has no tags.
- Rusanov or HLLC convective flux.
- Linear reconstruction from least-squares gradients with a Venkatakrishnan or Barth-Jespersen limiter, and a first-order fallback where the reconstruction would give negative density or pressure.
- Viscous fluxes from least-squares gradients with corrected face gradients, constant or Sutherland viscosity.
- Slip wall (Euler), or no-slip adiabatic / isothermal wall. Wall pressure from the limited face reconstruction.
- Characteristic (Riemann invariant) far-field boundary condition.
- SSP-RK2 / RK3; steady (local time step, including a viscous limit) or transient (global time step).
- CFL control with three modes: `off` (fixed CFL with a start-up ramp), `smart` (grows the CFL while steps are accepted, cuts it and rolls the step back on negative density or pressure, NaN, or a residual jump), `ser` (as `smart`, with the CFL also following the residual drop).
- Steady runs stop on a residual tolerance or when Cl, Cd and Cm vary by less than a set tolerance over a window.
- Angle of attack rotates the freestream; the mesh stays fixed.
- Lift, drag (pressure and friction), pitching moment about a chosen point (default quarter chord, positive nose-up), wall Cp and Cf, upper-surface separation point.
- Mean and half peak-to-peak amplitude of Cl, Cd and Cm over the last iterations.
- Plots: Mach number, pressure coefficient, vorticity, mesh near the body, wall Cp and Cf, and histories of residual, coefficients and CFL.
- Transient mode with live plots and optional frame saving.

---

## Test case: NACA 0012, laminar, M = 0.5, Re = 5000
- Mesh: `testcase_naca0012/naca0012.msh` from `gen_naca0012_gmsh.py` (12,302 cells: 8,396 triangles and 3,906 quads, boundary-layer quads, first cell height 6e-4, far field at 20 chords). The mesh passed `inspect_mesh.py` with no errors and 16 sliver warnings at the trailing edge.
- Reference: Swanson and Langer, NASA/TM-2016-219003 (finest grid): at alpha = 0 deg Cd 0.0556 (pressure 0.0228, friction 0.0328), separation at x/c 0.808; at alpha = 3 deg Cl 0.0513, Cd 0.0585 (pressure 0.0262, friction 0.0323), separation at x/c 0.464. The report gives no Cm.

| 14,000 iterations, mean of the last 4,000 (± half peak-to-peak) | alpha = 0 deg | alpha = 3 deg |
|---|---|---|
| Cl | 0.0002 ± 0.19 | 0.135 ± 0.26 |
| Cd | 0.0579 ± 0.0084 | 0.0632 ± 0.016 |
| Cm (about x = 0.25) | 0.0002 ± 0.089 | -0.014 ± 0.10 |

## Notes
- **The test case did not converge.** The density residual stalled (about 1e-2 at alpha = 0 deg, 0.3 to 0.45 at 3 deg) and Cl and Cm oscillate with large amplitude, so the solver does not reproduce the reference values. The mean Cd at alpha = 0 deg (0.0579) is about 4% above the reference, and friction drag (about 0.033 in the logs) is close to it, but the pressure drag and lift are not reliable. Separation points from the last iteration (0.667 and 0.405) are single snapshots of an oscillating flow.
- Likely causes, not confirmed: explicit Runge-Kutta stepping with local time steps is a weak solver for this case, the wake mesh is probably too coarse, and the very small fan cells at the trailing edge (area down to about 6e-8) slow the convergence. A real shedding instability cannot be ruled out.
- The CFL adaptation rejected no step in any run and sat at its cap, because the failure is a slow oscillation, not a positivity loss or residual jump. The rollback path itself has not been exercised. Lowering the CFL cap to 1 and raising the limiter constant to 100 at alpha = 0 deg did not change the behaviour. A test on stagnation of the residual, and implicit time stepping or implicit residual smoothing, are the next steps.
- Not tested: plotting, the Gmsh 4.x reader, the `mesh2d.py` text-file reader, the built-in O-grid path in this solver, isothermal wall, Sutherland viscosity, transient mode, and the Euler (`Re = Inf`) path. The cylinder at Re = 40 has not been re-run against the structured Version 5 solver, so do not treat this solver as validated.
- Runs take about 70 s per 1,000 iterations on the 12,300-cell mesh in GNU Octave 8.4 on one core.
- The angle of attack changes the freestream direction, not the mesh, so the built-in geometries are no longer rotated.
- Cm and the force coefficients use the reference length found from the wall nodes (the chord for an airfoil). Set `P.Lref` and `P.momentRef` to override.
- Gmsh files must be ASCII. Binary `.msh` files are rejected.
- Mesh quality is not checked inside the solver beyond zero or negative cell area and inconsistent orientation. Use `inspect_mesh.py` first.
- Grid convergence and validation are required before treating results as quantitatively reliable.

