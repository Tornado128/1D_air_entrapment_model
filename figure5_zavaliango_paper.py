## We generated figure 5 of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
## In this piece of code, a transient nonlinear partial differential equation has been solved using finite difference method to
# obtain the effect of punch velocity on the air pressure entrapped in the compact. This nonlinear equation (Eq. 7
# in Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012) is the result of incorporation of the continuum
# equation, ideal gas law and Darcy equation. The input parameters are initial gas pressure inside
# the powder, initial height of the compact, permeability and gas viscosity.
## Two punch velocities of 100 mm/s and 1000 mm/s are considered.

import numpy as np
import matplotlib.pyplot as plt
from material import *                                                      # constructor

Patm = 1                                                                    #atmospheric pressure (atm)
mu = 18.6*10**-6                                                            #viscosity of air at 25 C (kg/m.s)
N_L = 16                                                                    #number of numerical nodes in the compact

V = [100, 1000]                                                             #punch velocity (mm/s): we are assuming that both punches are moving at the same velocity
duration = [0.052, 0.0052]                                                  #duration of the simulation (second)

dt = 10**-6                                                                 #time step (second): time step should be as small as possible to avoid any numerical divergence
N_t = int(duration[0]/dt)+1                                                 #number of time steps for running the simulation. N_t will be shorter for high punch velocity
                                                                            # compared to low punch velocity. This is just an initialization

num_step = [int(duration[0]/dt)+1, int(duration[1]/dt)+1]                   #number of time steps for respectively low and high punch velocities

## Defining two dimensional lists for respectively height (H), relative density (D), porosity (phi) and permeability (K)
# vstime (the first column refers to the time and the second colum refers to punch velocity)
H = np.zeros((N_t, len(V)))
D = np.zeros((N_t, len(V)))
phi = np.zeros((N_t, len(V)))
K = np.zeros((N_t, len(V)))

## dz is the spatial increment in the vertical direction. As the compacts will be shrinking in size, the increment will
# be shrinking with time (m)
dz = np.zeros((N_t, len(V)))

H[0][0] = 8                                                                         #initial height of the compact under low velocity (mm)
H[0][1] = H[0][0]                                                                   #initial height of the compact under high velocity
D[0][0] = 0.33                                                                      #initial relative density of the compact under low velocity
D[0][1] = D[0][0]                                                                   #initial relative density of the compact under high velocity
phi[0][0] = 1 - D[0][0]                                                             #initial porosity of the compact under low velocity
phi[0][1] = phi[0][0]                                                               #initial porosity of the compact under high velocity

#print(r[0].a, r[0].b, r[0].c, r[0].d)
P = np.ones((N_t, N_L, len(V)))                                                     #defining and initilization of the air entrapment pressure to 1 atm for all the nodes at all the time steps

# looping over compacts (different punch velocities) and time steps
for k in range(len(V)):
    N_t = num_step[k]                                                               #The number of time steps is smaller for higher punch velocities
    for i in range(N_t):
        dz[i][k] = H[i][k]/N_L                                                      #increment size in z direction, which depends on time and punch velocities. It is getting smaller with time because the compacts are shrinking
        K [i][k] = r[0].a * np.exp(r[0].b*D[i][k])*abs(1-D[i][k]+r[0].c)**r[0].d
        if (i<N_t-1):
            H[i+1][k] = H[i][k] - V[k]*dt
            D[i+1][k] = D[0][k]*H[0][k]/H[i+1][k]
            phi[i+1][k] = 1 - D[i+1][k]
        for j in range(1, N_L):
            if (i<N_t-1):
                if (j<N_L-1):
                    A = P[i][j+1][k]
                    B = P[i][j][k]
                    C = P[i][j-1][k]


                    #finite difference scheme: Eq. (11) in Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                    P[i+1][j][k] = ( (1-phi[i+1][k])/phi[i+1][k] ) * ( (dt * K[i][k]*(10**5)*(A**2-2*B**2+C**2))/(2*mu*(1-phi[i][k])*dz[i][k]*dz[i][k]) + P[i][j][k]*phi[i][k]/(1-phi[i][k]) )
                else:
                    P[i+1][j][k] = P[i+1][j-1][k]

## fully entrapped air scenario
step = 0.01                                                                     #steps in porosity and relative density
D1 = np.arange(0.33, 0.92 + step, step).tolist()                              #discretizing from 0 to 0.8 with step size equals 'step'
size1 = len(D1)                                                                #size of the list
ratio1 = np.zeros(size1)
for i in range(size1):
    ratio1[i] = ( (1-np.min(D1))/(1-D1[i]) ) * D1[i]/np.min(D1)


#print(P.shape)
#selected_elements1 = P[:, 15, 0]
#selected_elements2 = P[:4894, 15, 1]
#print(selected_elements1)
#print(selected_elements2)
#print(D[:,0])
#print(D[:4894,1])

plt.semilogy(D[:num_step[0],0],P[:num_step[0], N_L-1, 0],'b--',D[:num_step[1],1],P[:num_step[1], N_L-1, 1],'r--',D1, ratio1,'g--')
plt.xlabel("relative density, D",fontsize=22)
plt.ylabel("P_max/P_atm",fontsize=22)
plt.ylim([10**0, 10**2])
plt.xlim([0.2, 1])
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(["punch velocity = 100 mm/s", "punch velocity = 1,000 mm/s","fully entrapped air"], fontsize=22)
plt.grid(True, which ="both")
plt.show()