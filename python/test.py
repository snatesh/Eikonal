from sFactors import *
from ngjQuad import *
from eikonal import *

import matplotlib.pyplot as plt

Np = 7
a = 0.5; b = a; c = a; d = a

Hab = sFactors(Np, a, b)
Habc = sFactors(Np, a, b, c)
Habcd = sFactors(Np, 0, 0, a, b, c, d) 
print(Hab)
print(Habc)
print(Habcd)

eikonal(0.5,0.5,0.5)

#Z0 = ngjQuad ( 5, 7, 0.5, 0.5, 0.5, 1e-14, 1e-14, 
#               0, 0, 0, 6)
#
#
#
#
#
#
#print(Z0)
#N = int(len(Z0)/3)
#
#
#
#
#
#
#plt.plot(Z0[0:N],Z0[N:2*N], 'o')
#plt.show()
               
