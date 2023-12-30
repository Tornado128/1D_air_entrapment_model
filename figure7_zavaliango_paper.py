## The goal in this piece of code is to generate figure 7 of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
## The goal is this piece of code is to obtain gas pressure profile inside a compact vs normalized distance from the
# center at 1000 mm/s at different final relative density. The maximum pressure in the case of symmetrical compression
# occurs at the center of the tablet

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from material import *                                                                      #constructor

Patm = 1                                                                                                        #atmospheric pressure (atm)
mu = 18.6*10**-6                                                                                                #viscosity of air at 25 C (kg/m.s)
N_L = 16                                                                                                        #number of numerical nodes in the compact

V = 1000                                                                                                        # punch velocity (mm/s)
dt = [8.131288E-7, 9.40128E-7, 0.000001]                                                                        # time step (s): The time steps are chosen in a way to produce the same number of time steps
duration = [0.00493, 0.0057, 0.006063]                                                                          # duration of the simulation (s)

## maximum pressure inside the compact, which occurs at the center
# of the compact in a symmetrical compaction
P_max = np.ones(len(dt))

## final height of the compacts after compression
final_height = np.ones(len(dt))
## final relative density of the compacts after compression
final_D = np.ones(len(dt))

N_t = int(duration[0]/dt[0]+1)                                                                                  #number of time steps for running the simulation. N_t will be shorter for high punch velocity
num_step = np.zeros(len(dt))
for i in range(len(duration)):
    num_step[i] = int(duration[i]/dt[i]+1)                                                                      #number of time steps for respectively low and high punch velocities

H_I = 8                                                                                                         #initial heigh of the compact (mm)
D_I = 0.23                                                                                                      #initial relative density
phi_I = 1 - D_I                                                                                                 #initial porosity
H = H_I * np.ones((N_t,len(dt)))                                                                                # height of the compact (mm)
height = H_I * np.ones((N_L,len(dt)))                                                                           # height of the compact (mm)
D = D_I * np.ones((N_t,len(dt)))                                                                                # relative density
phi = 1 - D                                                                                                     # porosity
K = np.ones((N_t,len(dt)))                                                                                      # permeability (mm^2)

## dz is the spatial increment in the vertical direction. As the compacts will be shrinking in size, the increment will
# be shrinking with time (m)
dz = np.zeros((N_t, len(dt)))

#defining and initilization of the air entrapment pressure (atm)
P = np.ones((N_t, N_L, len(dt)))

#initialization of the two-point normalized pressure profile in the compact
normalized_pressure = np.ones((N_L,len(dt)))

#initialization of the height of the compact
normalized_height = np.ones((N_L,len(dt)))


for k in range(len(dt)):
    print(k)
    N_t = int(num_step[k])                                                               #The number of time steps is smaller for higher punch velocities
    for i in range(int(N_t)):
        dz[i][k] = H[i][k]/N_L                                                      #increment size in z direction, which depends on time and punch velocities. It is getting smaller with time because the compacts are shrinking
        K [i][k] = r[0].a * np.exp(r[0].b*D[i][k])*abs(1-D[i][k]+r[0].c)**r[0].d    #permeability of Avicel-PH102
        if (i == int(N_t-1)):
            final_height[k] = H[i][k]                                               # final height of the compact (at the end of the compaction)
            final_D[k] = D[i][k]                                                    # final relative density (at the end of the compaction)
        if (i < int(N_t-1)):
            H[i+1][k] = H[i][k] - V*dt[k]
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

    P_max [k] = P[int(num_step[k]-1), N_L-1, k]                 # pressure center of the tablet for each final relative density

print(dz)
# obtaining the two-point normalized pressure
for k in range(len(dt)):
    N_t = int(num_step[k]) #The number of time steps is smaller for higher punch velocities
    for j in range(N_L):
        height[j][k] = ((N_L-1)-j)/(N_L-1)
        for i in range(int(N_t)):
            if (i==int(N_t)-1):
                normalized_pressure [j][k] = (P[N_t-1][j][k] - Patm)/(P_max[k]-Patm)

print(dz[:,2])
#changing the arrays to dataframe for the simplicity
dataframe_pressure = pd.DataFrame({'D_final=60%': normalized_pressure[:, 0], 'D_final=80%': normalized_pressure[:, 1], 'D_final=95%': normalized_pressure[:,2]})
print(dataframe_pressure)

dataframe_height = pd.DataFrame({'D_final=60%': height[:, 0], 'D_final=80%': height[:, 0], 'D_final=95%': height[:, 0]})
print(dataframe_height)

plt.plot(height[:,0],normalized_pressure[:,0],'b--',height[:,1],normalized_pressure[:,1],'r--',height[:,2],normalized_pressure[:,2],'g--',linewidth=3)
plt.xlabel("normalized distance from the center, z/(H/2)",fontsize=28)
plt.ylabel("(P - P_atm)/(P_max - P_atm)",fontsize=28)
plt.ylim([0, 1])
plt.xlim([0, 1])
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.grid(True, which ="both")
plt.legend(["D_final=60%", "D_final=80%","D_final=95%"], fontsize=28)
plt.show()


