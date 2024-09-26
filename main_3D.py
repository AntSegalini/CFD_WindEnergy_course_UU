from Fluid_3D import Fluid_windenergy
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.io import savemat,loadmat

# Parameters to set. 
Uinf, nu= 8, 10      # Parameters for the wind energy simulation.   
LX, NX= 1000, 100   # Parameters for the wind energy simulation.
LY, NY= 500, 64     # Parameters for the wind energy simulation.    
LZ, NZ= 500, 64     # Parameters for the wind energy simulation.

# turbine locations
X_T=[300]
Y_T=[0]
Z_T=[0]
R_T=[50]
C_T=[0.8]    
delta=10

t_max=LX/Uinf
iter_saving=100
start_from=0
CFL_max=0.2

fileName='simulation'

#########################################################################################

def main():

    sim = Fluid_windenergy(LX, NX, LY, NY, LZ, NZ, Uinf, nu, rho=1.225, X_T=X_T, Y_T=Y_T, Z_T=Z_T, R_T=R_T, C_T=C_T, delta=delta) # Create instance of the class.
    sim.compute_time_step(CFL_max)

    t,i=0,0
    if start_from>0:
        Q=loadmat(fileName+'_'+str(start_from)+'.mat')
        sim.set_uvwp(Q['u'],Q['v'],Q['w'],Q['p'])
        t=np.squeeze(Q['t'])
        i=start_from
        Q=0

    # Solve:
    t1=time.time()
    while (t<=t_max):
        i+=1
        t+=sim.dt
        
        sim.Momentum()
        sim.pressure_poisson_spectral()
        sim.corrector()
        sim.compute_time_step(CFL_max)
        
        print(i,' t=',np.round(t,1),' dt=',np.round(sim.dt,3),' MAX DIV=',np.round(np.max(np.abs(sim.DIVERGENCE(sim.u,sim.v,sim.w))),4))
        
        if np.mod(i,iter_saving)==0:
            savemat(fileName+'_'+str(i)+'.mat',{'u':sim.u,'v':sim.v,'w':sim.w,'p':sim.p,'t':t,'x_u':sim.x_u,'y_u':sim.y_u,'z_u':sim.z_u,'x_v':sim.x_v,'y_v':sim.y_v,'z_v':sim.z_v,'x_p':sim.x_p,'y_p':sim.y_p,'z_p':sim.z_p})

    savemat(fileName+'_'+str(i)+'.mat',{'u':sim.u,'v':sim.v,'w':sim.w,'p':sim.p,'t':t,'x_u':sim.x_u,'y_u':sim.y_u,'z_u':sim.z_u,'x_v':sim.x_v,'y_v':sim.y_v,'z_v':sim.z_v,'x_p':sim.x_p,'y_p':sim.y_p,'z_p':sim.z_p})

    print('elapsed time:',time.time()-t1)

    # # force volume integral (to calibrate a disk)
    # Thrust=-np.trapz(np.sum(np.sum(sim.alpha*sim.u**2,axis=0),axis=0))*sim.dx*sim.dy*sim.dz 
    # print('CT=',Thrust/(0.5*1.225*Uinf**2*np.pi*50**2)) 

    # part to plot something
    plt.subplot(221)    
    plt.contourf(sim.x_u,sim.y_u,sim.u[NZ//2],30,cmap=plt.get_cmap('jet'))
    plt.colorbar();plt.title('u');plt.xlabel('x [m]');plt.ylabel('y [m]')
    plt.subplot(222)
    plt.contourf(sim.x_v,sim.y_v,sim.v[NZ//2],30,cmap=plt.get_cmap('jet'))
    plt.colorbar();plt.title('v');plt.xlabel('x [m]');plt.ylabel('y [m]')
    plt.subplot(223)
    plt.contourf(sim.x_w,sim.y_w,sim.w[NZ//2],30,cmap=plt.get_cmap('jet'))
    plt.colorbar();plt.title('w');plt.xlabel('x [m]');plt.ylabel('y [m]')
    plt.subplot(224)
    plt.contourf(sim.x_p,sim.y_p,sim.p[NZ//2],30,cmap=plt.get_cmap('jet'))
    plt.colorbar();plt.title('p');plt.xlabel('x [m]');plt.ylabel('y [m]')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__": # Only run the functions defined in this code. 
    main()
