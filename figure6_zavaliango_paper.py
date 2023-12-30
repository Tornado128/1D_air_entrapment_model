## We generated figure 6 of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
## The goal is to obtain the maximum gas pressure vs the punch velocity.
## As the system and compaction are symmetrical, the maximum pressure occur in the middle.
## If the punch velocity is too high, an assymptotic behavior is observed, and the maximum entrapped gas pressure
# will approach to that of fully entrapped scenario

import numpy as np
import math
import matplotlib.pyplot as plt
from material import *                                                                                      #constructor: my database for permeability data
#from numpy import array

Patm = 1                                                                                                    #atmospheric pressure (atm)
mu = 18.6*10**-6                                                                                            #viscosity of air at 25 C (kg/m.s)
N_L = 16                                                                                                    #number of numerical nodes in the compact

V = [10, 100, 10**3, 10**4, 10**5, 10**6]                                                                   # punch velocity (mm/s)
dt = [10**-6, 10**-6,  10**-6, 10**-6, 10**-6, 10**-7]                                                      # time step (s)
duration = [0.5835, 0.05835, 0.005835, 0.0005835, 0.00005835, 0.000005835]                                  # duration of the simulation (s)

## maximum pressure inside the compact, which occurs at the center
# of the compact in a symmetrical compaction
P_max = np.ones(len(V))

N_t = int(duration[0]/dt[0]+1)                                                                              #number of time steps for running the simulation. N_t will be shorter for high punch velocity
num_step = np.zeros(len(V))
for i in range(len(duration)):
    num_step[i] = int(duration[i]/dt[i])+1                                                                  #number of time steps for respectively low and high punch velocities
H_I = 8                                                                                                     #initial heigh of the compact (mm)
D_I = 0.23                                                                                                  #initial relative density
phi_I = 1 - D_I                                                                                             #initial porosity
H = H_I * np.ones((N_t,len(V)))                                                                             # height of the compact (mm)
D = D_I * np.ones((N_t,len(V)))                                                                             # relative density
phi = 1 - D                                                                                                 # porosity
K = np.ones((N_t,len(V)))                                                                                   # permeability (mm^2)
## dz is the spatial increment in the vertical direction. As the compacts will be shrinking in size, the increment will
# be shrinking with time (m)
dz = np.zeros((N_t, len(V)))
#H[0] = 8                                                                                    # initial (t=0 second) height of the compact (mm)
#D[0] = 0.23                                                                                 # initial (t=0 second) relative density
#phi[0] = 1 - D[0]                                                                           # initial porosity

#print(r[0].a, r[0].b, r[0].c, r[0].d)
P = np.ones((N_t, N_L, len(V)))                                                     #defining and initilization of the air entrapment pressure to 1 atm for all the nodes at all the time steps

# looping over compacts (different punch velocities) and time steps
for k in range(len(V)):
    print(k)
    N_t = num_step[k]                                                               #The number of time steps is smaller for higher punch velocities
    for i in range(int(N_t)):
        dz[i][k] = H[i][k]/N_L                                                      #increment size in z direction, which depends on time and punch velocities. It is getting smaller with time because the compacts are shrinking
        K [i][k] = r[0].a * np.exp(r[0].b*D[i][k])*abs(1-D[i][k]+r[0].c)**r[0].d    #permeability of Avicel-PH102
        if (i<N_t-1):
            H[i+1][k] = H[i][k] - V[k]*dt[k]
            D[i+1][k] = D[0][k]*H[0][k]/H[i+1][k]
            phi[i+1][k] = 1 - D[i+1][k]
        for j in range(1, N_L):
            if (i<N_t-1):
                if (j<N_L-1):
                    A = P[i][j+1][k]
                    B = P[i][j][k]
                    C = P[i][j-1][k]


                    #finite difference scheme: Eq. (11) in Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                    P[i+1][j][k] = ( (1-phi[i+1][k])/phi[i+1][k] ) * ( (dt[k] * K[i][k]*(10**5)*(A**2-2*B**2+C**2))/(2*mu*(1-phi[i][k])*dz[i][k]*dz[i][k]) + P[i][j][k]*phi[i][k]/(1-phi[i][k]) )
                else:
                    P[i+1][j][k] = P[i+1][j-1][k]
    P_max [k] = P[int(num_step[k]-1), N_L-1, k]

print(P_max)
plt.loglog(V, P_max/Patm, 'ro',markersize=22)
plt.xlabel("punch velocity",fontsize=24)
plt.ylabel("P_max/P_atm",fontsize=24)
plt.xlim([10**-2, 10**6])
plt.ylim([1, 100])
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(True, which ="both")
plt.legend(["85% final relative density"], fontsize=22)
plt.show()