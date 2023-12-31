

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from material import *                                                                                          #constructor

Patm = 1                                                                                                        #atmospheric pressure (atm)
mu = 18.6*10**-6                                                                                                #viscosity of air at 25 C (kg/m.s)
N_L = 16                                                                                                        #number of numerical nodes in the compact

V = 1000                                                                                                        # punch velocity (mm/s)
d1 = 10                                                                                                         # punch diameter (mm)
L = d1 * 0.25                                                                                                   # length of punch die-gap (mm)
d2 = d1 + 40/1000                                                                                               # die diameter (mm)

area = 0.25*np.pi*(d1*d1)                                                                                       # punch area (mm^2)
dt = 10**-7                                                                                                     # time step (s)
duration = [0.00610, 0.0044]                                                                                    # duration of the simulation, which is different for the two materials (s)

N_t = int(duration[0]/dt+1)                                                                                     #number of time steps (longest) for running the simulation. N_t will be shorter for high punch velocity
num_step = np.zeros(len(duration))
for i in range(len(duration)):
    num_step[i] = int(duration[i]/dt+1)                                                                         #number of time steps for running for each height of the compact

H_I = [8, 8]                                                                                                  #initial heigh of the compacts (mm)
D_I = [r[0].e, r[1].e]                                                                                          #initial relative density for Avicel PH102 and Lactose 316
phi_I = [1-r[0].e, 1-r[1].e]                                                                                    #initial porosity
D = np.ones((len(duration),N_t))                                                                                #initial relative density
H = np.ones((len(duration),N_t))                                                                                #initialization of the height of the compact (mm)
K = np.ones((len(duration),N_t))                                                                                #initialization of the permeability (mm^2)
phi = np.zeros((len(duration),N_t))                                                                             #initialization of porosity


#defining and initilization of the air entrapment pressure (atm)
P = np.ones((len(duration),N_t, N_L))

for k in range(len(H_I)):
    H[k][0] = H_I[k]
    N_t = int(num_step[k])                                                                      #The number of time steps is smaller for higher punch velocities
    if (k==0):
        D[k][0]=D_I[k]
    else:
        D[k][0]=D_I[k]
    phi[k][0] = 1 - D[k][0]
    for i in range(int(N_t)):
        dz = H[k][i]/N_L                                                                        #increment size in z direction, which depends on time and punch velocities. It is getting smaller with time because the compacts are shrinking
        K[k][i] = r[k].a * np.exp(r[k].b*D[k][i])*abs(1-D[k][i]+r[k].c)**r[k].d                #permeability of Avicel-PH102
        if (i < int(N_t-1)):
            H[k][i+1] = H[k][i] - V*dt
            D[k][i+1] = D[k][0]*H[k][0]/H[k][i+1]
            phi[k][i+1] = 1 - D[k][i+1]
        for j in range(1, N_L):
            if (i<N_t-1):
                if (j<N_L-1):
                    A = P[k][i][j+1]
                    B = P[k][i][j]
                    C = P[k][i][j-1]

                    #finite difference scheme: Eq. (11) in Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                    P[k][i+1][j] = ( (1-phi[k][i+1])/phi[k][i+1] ) * ( (dt * K[k][i]*(10**5)*(A**2-2*B**2+C**2))/(2*mu*(1-phi[k][i])*dz*dz) + P[k][i][j]*phi[k][i]/(1-phi[k][i]) )
                    if (j == 1):
                        lambdai = d2 / d1  # die diameter tp punch diameter ratio
                        f = lambdai ** 4 - 1 - ((lambdai ** 2 - 1) ** 2) / np.log(lambdai)                                                  # see Eq. (18) of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                        alpha = np.pi * dz * f * d1 * d1 * d1 * d1 / (128 * K[k][i] * area * (1 - D[k][i]) * (L + H_I[k] - H[k][i]))           # see Eq. (19) of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                        P[k][i + 1][j - 1] = ((P[k][i + 1][j] ** 2 + alpha * Patm * Patm) / (1 + alpha)) ** 0.5                             # see Eq. (19) of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                else:
                    P[k][i+1][j] = P[k][i+1][j-1]

    #P_max [k] = P[int(num_step[k]-1), N_L-1,k]                 # pressure center of the tablet for each final relative density
    #print(H_I[k], P_max[k])

#print(D[1,:int(num_step[1])])
#dataframe_pressure = pd.DataFrame({'Avicel PH102': P[0,:,N_L-1], 'Lactose 316': P[1,:,N_L-1]})
#print(dataframe_pressure)

#dataframe_relative_density = pd.DataFrame({'Avicel PH102': D[0,:], 'Lactose 316': D[1,:]})
#print(dataframe_relative_density)

plt.semilogy(D[0,:int(num_step[0])], P[0,:int(num_step[0]),N_L-1], 'b-', D[1,:int(num_step[1])], P[1,:int(num_step[1]),N_L-1], 'r-',linewidth=3)
plt.xlabel("relative density",fontsize=24)
plt.ylabel("P_max/P_atm",fontsize=24)
plt.xlim([0, 1])
plt.ylim([1, 100])
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(True, which ="both")
plt.legend(["Avicel PH-102", "Lactose 316"], fontsize=28)
plt.show()
#dataframe_height = pd.DataFrame({'D_final=60%': height[:, 0], 'D_final=80%': height[:, 0], 'D_final=95%': height[:, 0]})
#print(dataframe_height)