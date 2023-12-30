## We generated figure 5 of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
## This piece of code estimates the permeability of two excipients vs relative density (D)
## The equations for permeability are given in Eq. (12) and Eq. (13) in "Zavaliangos, Antonios, et al.
## "Prediction of air entrapment in tableting: an approximate solution." Journal of pharmaceutical sciences 106.12 (2017):
# 3604-3612."

import numpy as np
import matplotlib.pyplot as plt
from material import *

D0 = 0.2                                                                            #initial relative density (=1-porosity) of the tablet
H0 = 16                                                                             #initial height of the tablet (mm)
V = 1000                                                                            #upper and lower punch velocity (mm/s)
step = 0.0001
T = 0.775*H0/(2*V)                                                                  #time of compression (s). 0.775 is just a number to make a maximum D of around 0.9
t = np.arange(0, T + step, step).tolist()                                           #discretizing from 0 to T with step size (s)

H = [H0 - i * 2*V for i in t]                                                       #Height of the compact vs time
D = [D0*H0/i for i in H]                                                            #relative density of the compact at different height

num_rows = len(r)
num_cols = len(D)
K = np.zeros((num_rows, num_cols))
for i in range(num_rows):
    for j in range(num_cols):
        K[i][j] = r[i].a * np.exp(r[i].b*D[j])*(1-D[j]+r[i].c)**r[i].d

plt.semilogy(D, K[0][:],'b--',D, K[1][:],'r--', D, K[2][:],'g--',D, K[3][:],'k--')
plt.xlabel("relative density, D",fontsize=22)
plt.ylabel("permeability,mm^2",fontsize=22)
plt.ylim([10**-10, 10**-4])
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(True, which ="both")
plt.legend(["%s"%r[0].name, "%s"%r[1].name, "%s"%r[2].name, "%s"%r[3].name], loc ="lower left", fontsize=22)
plt.show()



