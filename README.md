Finite difference code to simulate the flow around wind turbines in primitive variables (u,v,w,p) with a projection method to solve the incompressible Navier-Stokes equations.

The spanwise (Y) and vertical (Z) directions are periodic while a Dirichlet boundary condition is imposed in the velocity field at the inlet and a von Neumann boundary condition is imposed at the outlet. A constant eddy viscosity is adopted in the entire domain as turbulence closure.

To avoid the checkerboard problem, a staggered arrangement of the variables is adopted with a linear interpolation scheme in the convective term.
