# This constructor defines the properties of different materials including permeability and initial relative density
# The permeability equation is the same as Eq. (12) and Eq. (13) of Zavaliangos, Antonios, et al. "Prediction of air entrapment in tableting: an approximate solution." Journal of pharmaceutical sciences 106.12 (2017): 3604-3612.

 class material:
    def __init__(self, name, a, b, c, d, e):
        self.name = name
        self.a = a                                      #first parameter for permeability
        self.b = b                                      #second parameter for permeability
        self.c = c                                      #third parameter for permeability
        self.d = d                                      #fourth parameter for permeability
        self.e = e                                      #initial relative density

# K = a*exp(b*D)*(1-D-c)^d
# K is permeability and D is the relative density (=1-porosity): see equations 12 and 13 in
# Zavaliangos, Antonios, et al. "Prediction of air entrapment in tableting: an approximate solution."
# Journal of pharmaceutical sciences 106.12 (2017):
# 3604-3612."

n = 4               # number of materials that we have information about them
r =[0]*n            # defining a list with n elements of zero
r[0] = material("Avicel PH102 + 1% MgSt",3.7*10**-4,-10.49,-0.0312,2.088, 0.23)
r[1] = material("Lactose 316 + 1% MgSt",5.1*10**-2,-15.57,-0.0312,2.083, 0.435)
r[2] = material("VX-445 EG",5.6*10**-4,-10.16,-0.0999,2.144, 0.29)
r[3] = material("VX-121 EG",1.7*10**-4,-10.04,-0.12,1.739,0.23 )