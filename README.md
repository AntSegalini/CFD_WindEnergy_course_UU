Finite difference code to simulate the flow around wind turbines in primitive variables (u,v,w,p) with a projection method to solve the incompressible Navier-Stokes equations.

The spanwise (Y) and vertical (Z) directions are periodic while a Dirichlet boundary condition is imposed in the velocity field at the inlet and a von Neumann boundary condition is imposed at the outlet. A constant eddy viscosity is adopted in the entire domain as turbulence closure.

The spatial discretization is second-order accurate, while the temporal discretization is first-order accurate (Euler method). A pressure projection method is used to estimate at every timestep the pressure field that ensures the velocity field to be divergence-free. To avoid the checkerboard problem, a staggered arrangement of the variables is adopted with a linear interpolation scheme in the convective term.

Second-order discretizations are used to solve the momentum equation and to estimate the divergence. The pressure Poisson equation is solved as a spectral method in the periodic Y and Z directions so that only the discretization in the axial (X) direction is solved as a direct problem. Alternative formulations such as ADI techniques in Y and Z were not attempted.

The code can be run in python without any compilation. The reuired packages are numpy, scipy and matplotlib. The code is not parallelized so it runs with a single CPU. The user has to decide how long the simulation needs to be by setting a final time. It is possible to store the velocity and pressure field after some iterations (say every 100) and after the code reached the final time. Saving intermediate field is good to monitor the temporal evolution of the solution but these can be used to set a new initial condition and continue the simulation from that iteration in case the execution was interrupted or the solution did not reach a steady-state yet.

The user can set the parameters of the simulation at the beginning of the main_3d.py file as
 
Uinf, nu= 8, 10      # free-stream velocity [m/s], eddy viscosity [m2/s]

LX, NX= 1000, 100    # X domain length [m], number of grid points in X (x=0 is the domain start)

LY, NY= 500, 64      # Y domain length [m], number of grid points in Y (periodic direction) (y=0 is the domain center)

LZ, NZ= 500, 64      # Z domain length [m], number of grid points in Z (periodic direction) (z=0 is the domain center)

X_T=[300]            # X location of the turbine [m] (for two turbines X_T=[300,500])

Y_T=[0]              # Y location of the turbine [m] (for two turbines Y_T=[0,10])

Z_T=[0]              # Z location of the turbine [m] (for two turbines Z_T=[0,0])

R_T=[50]             # Radius of the turbine [m] (for two turbines R_T=[50,50])

C_T=[0.8]            # Thrust coefficient of the turbine (for two turbines C_T=[0.8,0.7])

delta=10             # coefficient used to smear in the axial direction the turbine thrust force as a Gaussian force distribution with standard deviation delta (it should be comparable to Lx/Nx)

t_max=LX/Uinf        # Final simulation time [s] (it can be a number)

iter_saving=100      # number of iteration after which the actual velocity and pressure fields are stored

start_from=0         # starting iteration. If it is zero, the field is initialized as the free-stream velocity everywhere, if non zero, it loads the field stored previously at the iteration start_from

CFL_max=0.2          # Maximum allowed CFL number

fileName='simulation'  # Name of the stored fields


When the final time is reached, the code stores the field in a Matlab file (this is the same for the intermediate steps as well) and it plots some hub-height planes to evaluate the solution quality. Each Matlab file contains the following variables

   Attr Name        Size                     Bytes  Class
   
   ==== ====        ====                     =====  ===== 
   
        p          64x64x99                3244032  double
        
        t           1x1                          8  double
        
        u          64x64x100               3276800  double
        
        v          64x64x101               3309568  double
        
        w          64x64x101               3309568  double
        
        x_p         1x99                       792  double
        
        x_u         1x100                      800  double
        
        x_v         1x101                      808  double
        
        y_p         1x64                       512  double
        
        y_u         1x64                       512  double
        
        y_v         1x64                       512  double
        
        z_p         1x64                       512  double
        
        z_u         1x64                       512  double
        
        z_v         1x64                       512  double

The 3D arrays for the velocity components (u,v,w) in the (x,y,z) directions and the pressure are arranged so that the first index is for the Z dimension, the second index for the Y dimension and the third index is for the X dimension. Since the variables are defined in a staggered grid, the grid coordinates are also provided. For instance, for the axial velocity u, the X dimension is discretized by x_u, the Y dimension is discretized by y_u and so forth. The user is reminded that u,v,w and p might have different size since they are defined at different grid points

