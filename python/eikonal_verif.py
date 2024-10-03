import numpy as np
from jPoly import *
from pMat import *
import matplotlib.pyplot as plt

Z = np.loadtxt("triquadLeg_15_24.txt"); Nq = int(len(Z) / 3.0)
Xq = Z[0:Nq]; Yq = Z[Nq:2*Nq]; Wq = Z[2*Nq:3*Nq]
#####################################################################

n = 15; a = b = c = 0.5; nthreads = 6
xx, yy = np.meshgrid(np.linspace(0,1,30), np.linspace(0,1,30))
X = xx[yy[:]<1-xx[:]]; Y = yy[yy[:]<1-xx[:]]
Nx = len(X); 
Np = int(0.5*(n+1)*(n+2))

def u(x1, x2):
  dists = np.zeros_like(x1)
  for j in range(np.size(x1,0)):
    x = np.array([x1[j],x2[j]])
    z = np.array([(x1[j]-x2[j]+1.0)/2.0, (x2[j]-x1[j]+1.0)/2.0])
    d = np.linalg.norm(x-z)
    dists[j] = np.min([x1[j],x2[j],d])
  return dists

def f(x, y):
  return np.sin(x**2 + y**2)

def dxf(x,y):
  return 2*x*np.cos(x**2 + y**2)

def dyf(x,y):
  return 2*y*np.cos(x**2 + y**2)


abcPoly = jPoly(Nx, n, a, b, c, nthreads)
Vabc = computeV(abcPoly, Np, X, Y)
a1bc1Poly = jPoly(Nx, n, a+1, b+1, c+1, nthreads)
Va1bc1 = computeV(a1bc1Poly, Np, X, Y)

cf, Vabcq = computeCoeffs(f, n, a, b, c, Xq, Yq, Wq, nthreads)
print(np.size(cf))
print(np.linalg.norm(Wq.dot(Vabcq.dot(cf) - f(Xq,Yq))))



#####################################################################



cu, _ = computeCoeffs(u, n, a, b, c, Xq, Yq, Wq, nthreads)



fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, Y, u(X,Y))
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, Y, V.dot(cu))
plt.show()

