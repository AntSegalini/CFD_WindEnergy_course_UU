import numpy as np
from scipy.linalg import solve_banded
# from concurrent.futures import ProcessPoolExecutor

def process_i(i, kz, ky, D, I, b_hat, jmin=0):
    result=np.zeros((len(ky),b_hat.shape[-1]),dtype=np.complex128)
    for j in range(jmin,len(ky)):
        result[j] = solve_banded((1, 1), D - (kz**2 + ky[j]**2) * I, b_hat[j])
    return (i, result)

class Fluid_windenergy:
    # staggered grid version
    def __init__(self, LX, NX, LY, NY, LZ, NZ, Uinf, nu, rho=1.225, X_T=[], Y_T=[], Z_T=[], R_T=[], C_T=[], delta=0.1):
        self.NX = NX
        self.NY = NY
        self.NZ = NZ

        self.LX = LX
        self.LY = LY
        self.LZ = LZ

        # Initialize arrays. 
        self.u = np.zeros((NZ, NY, NX)) + Uinf    # dy/2, dz/2 (no ghost cells)
        self.v = np.zeros((NZ, NY, NX+1))         # dx/2, dz/2 (2 ghost cells)
        self.w = np.zeros((NZ, NY, NX+1))         # dx/2, dy/2 (2 ghost cells)    
        self.p = np.zeros((NZ, NY, NX-1))         # dx/2, dy/2, dz/2 (no ghost cells)
        self.alpha= np.zeros_like(self.u)
        self.Uinf = Uinf

        # Define steps 
        x= np.linspace(0, LX, NX)
        y= np.linspace(0, LY, NY, endpoint=False)-LY//2
        z= np.linspace(0, LZ, NZ, endpoint=False)-LZ//2

        self.dx = LX / (NX - 1)
        self.dy = LY / NY   # periodic direction
        self.dz = LZ / NZ   # periodic direction
        
        self.x_u=x
        self.y_u=y+self.dy/2
        self.z_u=z+self.dz/2
        
        self.x_v=np.linspace(-self.dx/2, LX+self.dx/2, NX+1)
        self.y_v=y
        self.z_v=z+self.dz/2

        self.x_w=self.x_v
        self.y_w=y+self.dy/2
        self.z_w=z

        self.x_p=self.x_v[1:-1]
        self.y_p=y+self.dy/2
        self.z_p=z+self.dz/2

        self.nu = nu
        self.compute_time_step(0.7)

        self.D=np.ones((3,NX-1))/self.dx**2
        self.D[1]*=-2
        self.D[1,[0,-1]]/=2         # zero pressure derivative at outlet/inlet
        
        # Boundary conditions for pressure for zeroth wavenumber
        self.D0=self.D.copy()
        self.D0[1,0]*=3             # zero pressure at inlet
        
        self.I=np.zeros((3,NX-1))
        self.I[1]=1
        
        # wavenumbers to be used in the Poisson solver
        self.ky=np.fft.rfftfreq(NY)*(2*np.pi/self.dy)
        self.kz=np.fft.fftfreq(NZ)*(2*np.pi/self.dz)

        Z,Y=np.meshgrid(self.z_u,self.y_u,indexing='ij') # the shift is added to have the correct location of the grid points for the X forcing
        xp=np.zeros_like(x)
        rp=np.zeros_like(Y)
        for i in range(len(X_T)):
            a=(1-np.sqrt(1-C_T[i]))/2   # result from actuator disk theory (to be calibrated)

            rp[:]=((Y-Y_T[i])/R_T[i])**2 + ((Z-Z_T[i])/R_T[i])**2
            x_pos=np.where(np.abs(x-X_T[i])<3*delta)[0]
            xp[x_pos]=np.exp(-((x[x_pos]-X_T[i])/delta)**2/2)   # Gaussian streamwise smearing
            xp/=np.trapz(xp)*self.dx                            # normalization of the Gaussian smearing
            r_pos=np.where(rp<1.3225)
            r_pos2=np.where(rp>=1.3225)
            rp[r_pos]=1+np.tanh((1-rp[r_pos])/0.05)
            rp[r_pos2]=0
            factor=np.trapz(np.trapz(rp,axis=0),axis=0)*self.dy*self.dz/(np.pi*R_T[i]**2)
            rp[r_pos]/=factor
            for j in x_pos:
                self.alpha[:,:,j]+= -rp * xp[j] *0.5*rho*C_T[i]/(1-a)**2
            xp[x_pos]=0

    def set_uvwp(self,u,v,w,p):
        self.u[:]=u.copy()
        self.v[:]=v.copy()
        self.w[:]=w.copy()
        self.p[:]=p.copy()

    def X_derivative(self, p, scheme='grid'):
        if scheme=='grid':
            # derivative at the grid point (second-order)
            return (p[:,:, 2:] - p[:,:, :-2]) / (2 * self.dx)
        elif scheme=='middle':
            # derivative between two grid points (second-order)
            return (p[:,:, 1:] - p[:,:, :-1]) / self.dx 
        
    def Y_derivative(self, p, scheme='grid'):
        if scheme=='grid':
            # derivative at the grid point (second-order)
            return (np.roll(p,-1,axis=1) - np.roll(p,1,axis=1)) / (2 * self.dy) 
        elif scheme=='middle':
            # derivative between two grid points (second-order)
            return (np.roll(p,-1,axis=1) - p) / self.dy
            
    def Z_derivative(self, p, scheme='grid'):
        if scheme=='grid':
            # derivative at the grid point (second-order)
            return (np.roll(p,-1,axis=0) - np.roll(p,1,axis=0)) / (2 * self.dz) 
        elif scheme=='middle':
            # derivative between two grid points (second-order)
            return (np.roll(p,-1,axis=0) - p) / self.dz    

    def LAPLACIAN(self, p):
        # Laplacian at the internal points (x)
        q=p[:,:, 1:-1]
        return  ((p[:,:, 2:] - 2 * q + p[:,:, :-2]) / self.dx**2 +
         (np.roll(q,-1,axis=1) - 2 * q + np.roll(q,1,axis=1) ) / self.dy**2 +
         (np.roll(q,-1,axis=0) - 2 * q + np.roll(q,1,axis=0) ) / self.dz**2 )
    
    def DIVERGENCE(self,u,v,w):
        return self.X_derivative(u,'middle') + self.Y_derivative(v[:,:, 1:-1],'middle') + self.Z_derivative(w[:,:, 1:-1],'middle')

    def INTERPOLATE_u(self,p,variable='v'):
        q=p[:,:, 1:-2] + p[:,:, 2:-1]       # averaging in x direction (both v and w)
        if variable=='v':
            q += np.roll(q,-1,axis=1)       # averaging in y direction
        elif variable=='w':
            q += np.roll(q,-1,axis=0)       # averaging in z direction
        q/=4
        return q
    
    def INTERPOLATE_v(self,p,variable='u'):
        if variable=='u':
            q=p[:,:, 1:]+p[:,:, :-1]        # averaging in x direction
        elif variable=='w':
            q=p[:,:, 1:-1]+np.roll(p[:,:, 1:-1],-1,axis=0)      # averaging in z direction
        q+=np.roll(q,1,axis=1)              # averaging in y direction (both u and w)
        q/=4
        return q
    
    def INTERPOLATE_w(self,p,variable='u'):
        if variable=='u':
            q=p[:,:, 1:]+p[:,:, :-1]        # averaging in x direction
        elif variable=='v':
            q=p[:,:, 1:-1]+np.roll(p[:,:, 1:-1],-1,axis=1)    # averaging in y direction
        q+=np.roll(q,1,axis=0)              # averaging in z direction (both u and v)
        q/=4
        return q
##############################################################################

    def Momentum (self):
        # Initialize arrays 
        u, v, w, dt, nu= self.u, self.v, self.w, self.dt, self.nu
        un = u.copy()
        vn = v.copy()
        wn = w.copy()

        # This excludes the inlet/outlet boundaries. 
        u[:,:, 1:-1] = un[:,:, 1:-1] + dt*(
                        - un[:,:, 1:-1]              * self.X_derivative(un)
                        - self.INTERPOLATE_u(vn,'v') * self.Y_derivative(un[:,:, 1:-1])
                        - self.INTERPOLATE_u(wn,'w') * self.Z_derivative(un[:,:, 1:-1])
                        + nu * self.LAPLACIAN(un)
                        + self.alpha[:,:,1:-1]*un[:,:,1:-1]**2)
        
        v[:,:, 1:-1] = vn[:,:, 1:-1] + dt*(
                        - self.INTERPOLATE_v(un,'u') * self.X_derivative(vn)
                        - vn[:,:, 1:-1]              * self.Y_derivative(vn[:,:, 1:-1])
                        - self.INTERPOLATE_v(wn,'w') * self.Z_derivative(vn[:,:, 1:-1])
                        + nu * self.LAPLACIAN(vn) )
                
        w[:,:, 1:-1] = wn[:,:, 1:-1] + dt*(
                        - self.INTERPOLATE_w(un,'u') * self.X_derivative(wn)
                        - self.INTERPOLATE_w(vn,'v') * self.Y_derivative(wn[:,:, 1:-1])
                        - wn[:,:, 1:-1]              * self.Z_derivative(wn[:,:, 1:-1])
                        + nu * self.LAPLACIAN(wn) )
        
        # Inlet, Dirichlet BC.
        u[:,:, 0] = self.Uinf
        v[:,:, 0] = -v[:,:, 1] # ghost cell (v=0)
        w[:,:, 0] = -w[:,:, 1] # ghost cell (w=0)
        # Outlet, Neumann. du/dx=0.
        u[:,:, -1] = (4*u[:,:, -2]-u[:,:, -3])/3
        v[:,:, -1] = v[:,:, -2]
        w[:,:, -1] = w[:,:, -2]
        # periodic BC on top and bottom (no action required)
        
    def pressure_poisson_spectral(self):
        u, v, w, dt = self.u, self.v, self.w, self.dt
        
        b_hat=np.fft.rfft2( 1 / dt * self.DIVERGENCE(u, v, w) ,axes=(0,1))
        
        b_hat[0,0] = solve_banded((1, 1), self.D0, b_hat[0, 0])
        b_hat[0,1:] = process_i(0, 0, self.ky, self.D, self.I, b_hat[0],1)[1][1:]

        for i in range(1,len(self.kz)):
            b_hat[i] = process_i(i, self.kz[i], self.ky, self.D, self.I, b_hat[i])[1]

        self.p[:]=np.fft.irfft2(b_hat,axes=(0,1))
        
    def corrector(self):
        u, v, w= self.u, self.v, self.w
        # Corrector step. 
        u[:,:, 1:-1] -= self.dt * self.X_derivative(self.p,'middle')
        v[:,:, 1:-1] -= self.dt * self.Y_derivative(np.roll(self.p,1,axis=1),'middle')
        w[:,:, 1:-1] -= self.dt * self.Z_derivative(np.roll(self.p,1,axis=0),'middle')
        # Inlet, Dirichlet BC.
        u[:,:, 0] = self.Uinf
        v[:,:, 0] = -v[:,:, 1] # ghost cell (v=0)
        w[:,:, 0] = -w[:,:, 1] # ghost cell (w=0)
        # Outlet, Neumann BC. du/dx=0.
        u[:,:, -1] = (4*u[:,:, -2]-u[:,:, -3])/3
        v[:,:, -1] = v[:,:, -2]
        w[:,:, -1] = w[:,:, -2]
        # periodic BC on top and bottom (no action required)

    def compute_time_step(self, cfl_number=0.1):
        """
        Compute the maximum allowable time step for a 3D incompressible Navier-Stokes simulation.
        I could use this function to automatically adjust the timestep based on the current velocity field. 
        
        Parameters:
            cfl_number (float): CFL number (typically < 1 for stability).
            
        Returns:
            float: Maximum allowable time step.
        """
        u, v, w, dx, dy, dz, nu= self.u, self.v, self.w, self.dx, self.dy, self.dz, self.nu
        u_max = np.max(np.abs(u))
        v_max = np.max(np.abs(v))
        w_max = np.max(np.abs(w))
        
        dt_conv_x = dx / u_max if u_max != 0 else np.inf
        dt_conv_y = dy / v_max if v_max != 0 else np.inf
        dt_conv_z = dz / w_max if w_max != 0 else np.inf
        dt_diff_x = dx**2 / (2 * nu)
        dt_diff_y = dy**2 / (2 * nu)
        dt_diff_z = dz**2 / (2 * nu)
        
        self.dt = cfl_number * min(dt_conv_x, dt_conv_y, dt_conv_z, dt_diff_x, dt_diff_y, dt_diff_z)
        