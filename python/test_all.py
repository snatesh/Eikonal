from pMat import *
from sFactors import *
from jPoly import *

# jacobi poly params
a = 0.5; b = 0.5; c = 0.5;
Z = np.loadtxt("triquadleg_16_27.txt")
N = int(len(Z) / 3.0)
print(N)
