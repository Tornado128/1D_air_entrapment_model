## This piece of code generates figure 8 in Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
## This piece of code explores the effect of punch-die gap on the entrapped air pressure inside a compact

import matplotlib.pyplot as plt
import numpy as np
from material import *                                                              # constructor

Patm = 1                                                                            #atmospheric pressure (atm)
d1 = 8.0                                                                            #punch diameter (mm)
L = d1 * 0.21                                                                       #Initial vertical length of of punch-die gap
mu = 18.6*10**-6                                                                    #viscosity of air at 25 C (kg/m.s)
N_L = 16                                                                            #number of numerical nodes in the compact

dt = 3*10**-8                                                                       #time step (second): time step should be as small as possible to avoid any numerical divergence
V = 1000                                                                            #punch velocity (mm/s): we are assuming that both punches are moving at the same velocity
duration = 0.00596                                                                  #duration of the simulation (second)
d2 = np.arange(d1 + 0.0000001, d1 + 0.150, 0.005).tolist()                          #die diameter (mm)
area = np.pi*d1*d1/4                                                                #cross sectional area of the tablet

# gap between the punch and die (mm)
gap = np.zeros(len(d2))


N_t = int(duration/dt)+1                                                            #number of time steps for running the simulation. N_t will be shorter for high punch velocity
                                                                                    # compared to low punch velocity. This is just an initialization
##initialization of the height of the compact
P = np.ones((len(d2),N_t,N_L ))

##Defining one-dimensional lists for respectively height (H), relative density (D), porosity (phi) and permeability (K) vstime
H_I = 8                                                                                 #initial height of the compact under low velocity (mm)
D_I = r[0].e                                                                              #initial relative density of the compact under low velocity

H = H_I * np.ones(N_t)
D = D_I * np.ones(N_t)
phi = 1 - D
dz = np.zeros(N_t)
K = np.zeros(N_t)                                                                   #itinializing permeability (mm^2)

## dz is the spatial increment in the vertical direction. As the compacts will be shrinking in size, the increment will
# be shrinking with time (m)
dz = np.zeros(N_t)

## maximum pressure inside the compact, which occurs at the center
# of the compact in a symmetrical compaction
P_max = np.ones(len(d2))

# looping over compacts (different punch velocities) and time steps
for k in range(len(d2)):
    gap[k] = 1000 * (d2[k] - d1)                                                            #gap between the punch and die (micron)
    for i in range(N_t):
        dz[i] = H[i]/N_L                                                                    #increment size in z direction, which depends on time and punch velocities. It is getting smaller with time because the compacts are shrinking
        K [i] = r[0].a * np.exp(r[0].b*D[i])*abs(1-D[i]+r[0].c)**r[0].d                     #permeability of Avicel-PH102
        if (i<N_t-1):
            H[i+1] = H[i] - V * dt
            D[i+1] = D[0]*H[0]/H[i+1]
            phi[i+1] = 1 - D[i+1]
        for j in range(1, N_L):
            if (i<N_t-1):
                if (j<N_L-1):
                    A = P[k][i][j+1]
                    B = P[k][i][j]
                    C = P[k][i][j-1]

                    #finite difference scheme: Eq. (11) in Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                    P[k][i+1][j] = ( (1-phi[i+1])/phi[i+1] ) * ( (dt * K[i]*(10**5)*(A**2-2*B**2+C**2))/(2*mu*(1-phi[i])*dz[i]*dz[i]) + P[k][i][j]*phi[i]/(1-phi[i]) )
                    if (j == 1):
                        lambdai = d2[k]/d1                                                          #die diameter tp punch diameter ratio
                        f = lambdai**4 - 1 - ((lambdai**2-1)**2)/np.log(lambdai)                    #see Eq. (18) of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                        alpha = np.pi*dz[i]*f*d1*d1*d1*d1/(128*K[i]*area*(1-D[i])*(L+H_I-H[i]))     #see Eq. (19) of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                        P[k][i+1][j-1] = ((P[k][i+1][j]**2 + alpha*Patm*Patm)/(1+alpha))**0.5       #see Eq. (19) of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
                else:
                    P[k][i+1][j] = P[k][i+1][j-1]
    P_max [k] = P[k, N_t-1, N_L-1]
    print(k, P_max[k], D[N_t-1], gap[k])

plt.plot(gap, P_max/Patm, 'ro',markersize=18)
plt.xlabel("punch-die gap, micron",fontsize=24)
plt.ylabel("P_max/P_atm",fontsize=24)
plt.xlim([0, np.max(gap)])
plt.ylim([0, 35])
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(True, which ="both")
plt.legend(["90% final relative density"], fontsize=22)
plt.show()