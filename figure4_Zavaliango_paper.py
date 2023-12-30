## We generated figure 4 of Zavaliangos et al, J. Pharmaceutical Sciences, 106 (2017) 3604-36012
## This piece of code estimates the pressure of FULLY entrapped air as a fuction of compact relative density (D) for
# 3 different levels of initial relative density (D). It is based on Boyle's law because air can not leave the compact
## The equations for permeability are given in Eq. (12) and Eq. (13) in "Zavaliangos, Antonios, et al.
# "Prediction of air entrapment in tableting: an approximate solution." Journal of pharmaceutical sciences 106.12 (2017):
# 3604-3612."

import matplotlib.pyplot as plt
from material import *
import numpy as np

P0 = 1                                                                         #atmospheric pressure (atm)
step = 0.05                                                                    #steps in porosity and relative density


phi1 = np.arange(0.05, 0.8 + step, step).tolist()                              #discretizing from 0 to 0.8 with step size equals 'step'
D1 = [1-i for i in phi1]                                                       #relative density
size1 = len(D1)                                                                #size of the list

phi2 = np.arange(0.05, 0.7 + step, step).tolist()                              #discretizing from 0 to 0.7 with step size equals 'step'
D2 = [1-i for i in phi2]                                                       #relative density
size2 = len(D2)                                                                #size of the list

phi3 = np.arange(0.05, 0.6 + step, step).tolist()                              #discretizing from 0 to 0.6 with step size equals 'step'
D3 = [1-i for i in phi3]                                                       #relative density
size3 = len(D3)

ratio1 = np.zeros(size1)
for i in range(size1):
    ratio1[i] = ( (1-np.min(D1))/(1-D1[i]) ) * D1[i]/np.min(D1)
ratio2 = np.zeros(size2)
for i in range(size2):
    ratio2[i] = ( (1-np.min(D2))/(1-D2[i]) ) * D2[i]/np.min(D2)
ratio3 = np.zeros(size3)
for i in range(size3):
    ratio3[i] = ( (1-np.min(D3))/(1-D3[i]) ) * D3[i]/np.min(D3)

plt.semilogy(D1, ratio1,'b--',D2, ratio2,'r--', D3, ratio3, 'g--')
plt.xlabel("relative density, D",fontsize=22)
plt.ylabel("P_{max}/P_{atm}",fontsize=22)
plt.ylim([10**0, 10**2])
plt.xlim([0,1])
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(True, which ="both")
plt.legend(["D_0=0.2","D_1=0.3","D_2=0.4"], fontsize=22)
plt.show()



print(ratio1)



