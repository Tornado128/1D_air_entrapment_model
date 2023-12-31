## This piece of code generates figure 9 in Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
## This piece of code explores the effect of the initial die filling on the air pressure in the middle of the compact

import numpy as np
import matplotlib.pyplot as plt
from material import *                                                                                          #constructor

Patm = 1                                                                                                        #atmospheric pressure (atm)
mu = 18.6*10**-6                                                                                                #viscosity of air at 25 C (kg/m.s)
N_L = 16                                                                                                        #number of numerical nodes in the compact


V = 1000                                                                                                        # punch velocity (mm/s)
d1 = 17.5                                                                                                       # punch diameter (mm)
L = d1 * 0.21                                                                                                   # length of punch die-gap (mm)
area = 0.25*np.pi*(d1*d1)                                                                                       # punch area (mm^2)
d2 = d1 + 53/1000                                                                                               # die diameter (mm)
dt = 10**-7                                                                                                     # time step (s): The time steps are chosen in a way to produce the same number of time steps
duration = [0.011911, 0.010422, 0.008933, 0.008189, 0.007445, 0.0067, 0.005956]                                # duration of the simulation (s)

N_t = int(duration[0]/dt+1)                                                                                     #number of time steps (longest) for running the simulation. N_t will be shorter for high punch velocity
num_step = np.zeros(len(duration))
for i in range(len(duration)):
    num_step[i] = int(duration[i]/dt+1)                                                                         #number of time steps for running for each height of the compact

H_I = [16, 14, 12, 11, 10, 9, 8]                                                                                #initial heigh of the compacts (mm)
D_I = r[0].e * np.ones((N_t,len(duration)))                                                                     #initial relative density
phi_I = 1 - D_I                                                                                                 #initial porosity
H = np.ones((N_t,len(duration)))                                                                                # initialization of the height of the compact (mm)
D = D_I * np.ones((N_t,len(duration)))                                                                          # initialization of the relative density
phi = 1 - D                                                                                                     # initialization of the porosity
K = np.ones((N_t,len(duration)))                                                                                # initialization of the permeability (mm^2)

#defining and initilization of the air entrapment pressure (atm)
P = np.ones((N_t, N_L, len(duration)))

## maximum pressure inside each compact, which occurs at the center
# of the compact in a symmetrical compaction
P_max = np.ones(len(H_I))

for k in range(len(H_I)):
    H[0][k] = H_I[k]
    N_t = int(num_step[k])                                                                      #The number of time steps is smaller for higher punch velocities
    for i in range(int(N_t)):
        dz = H[i][k]/N_L                                                                        #increment size in z direction, which depends on time and punch velocities. It is getting smaller with time because the compacts are shrinking
        K [i][k] = r[0].a * np.exp(r[0].b*D[i][k])*abs(1-D[i][k]+r[0].c)**r[0].d                #permeability of Avicel-PH102
        if (i < int(N_t-1)):
            H[i+1][k] = H[i][k] - V*dt
            D[i+1][k] = D[0][k]*H[0][k]/H[i+1][k]
            phi[i+1][k] = 1 - D[i+1][k]
        for j in range(1, N_L):
            if (i<N_t-1):
                if (j<N_L-1):
                    A = P[i][j+1][k]
                    B = P[i][j][k]
                    C = P[i][j-1][k]

                    #finite difference scheme: Eq. (11) in Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                    P[i+1][j][k] = ( (1-phi[i+1][k])/phi[i+1][k] ) * ( (dt * K[i][k]*(10**5)*(A**2-2*B**2+C**2))/(2*mu*(1-phi[i][k])*dz*dz) + P[i][j][k]*phi[i][k]/(1-phi[i][k]) )
                    if (j == 1):
                        lambdai = d2 / d1  # die diameter tp punch diameter ratio
                        f = lambdai ** 4 - 1 - ((lambdai ** 2 - 1) ** 2) / np.log(lambdai)                                                  # see Eq. (18) of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                        alpha = np.pi * dz * f * d1 * d1 * d1 * d1 / (128 * K[i][k] * area * (1 - D[i][k]) * (L + H_I[k] - H[i][k]))           # see Eq. (19) of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                        P[i + 1][j - 1][k] = ((P[i + 1][j][k] ** 2 + alpha * Patm * Patm) / (1 + alpha)) ** 0.5                             # see Eq. (19) of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                else:
                    P[i+1][j][k] = P[i+1][j-1][k]

    P_max [k] = P[int(num_step[k]-1), N_L-1,k]                 # pressure center of the tablet for each final relative density
    print(H_I[k], P_max[k])

plt.plot(H_I, P_max/Patm, 'ro',markersize=18)
plt.xlabel("initial die filling height, micron (H/2)",fontsize=24)
plt.ylabel("P_max/P_atm",fontsize=24)
plt.xlim([0, np.max(H_I)])
plt.xlim([16, 24])
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(True, which ="both")
plt.legend(["90% final relative density"], fontsize=22)
plt.show()