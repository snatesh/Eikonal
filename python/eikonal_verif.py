import numpy as np
from jPoly import *
from pMat import *
import matplotlib.pyplot as plt

#####################################################################

Z = np.loadtxt("triquadleg_n22_m34_N276.txt"); Nq = int(len(Z) / 3.0)
Xq = Z[0:Nq]; Yq = Z[Nq:2*Nq]; Wq = 2.0 * Z[2*Nq:3*Nq]

def u(x1, x2):
  dists = np.zeros_like(x1)
  for j in range(np.size(x1,0)):
    x = np.array([x1[j],x2[j]])
    z = np.array([(x1[j]-x2[j]+1.0)/2.0, (x2[j]-x1[j]+1.0)/2.0])
    d = np.linalg.norm(x-z)
    dists[j] = np.min([x1[j],x2[j],d])
  return dists
#####################################################################

n = 10
a = b = c = 0.5
nthreads = 6

cu, _ = computeCoeffs(u, n, a, b, c, Xq, Yq, Wq, nthreads)

xx, yy = np.meshgrid(np.linspace(0,1,30), np.linspace(0,1,30))
X = xx[yy[:]<1-xx[:]]; Y = yy[yy[:]<1-xx[:]]
Nx = len(X); N = int(0.5*(n+1)*(n+2))
poly = jPoly(Nx, n, a, b, c, nthreads)
V = computeV(poly, N, X, Y)

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, Y, u(X,Y))
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, Y, V.dot(cu))
plt.show()

