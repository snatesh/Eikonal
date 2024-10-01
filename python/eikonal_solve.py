from numpy import *
from eikonal import *
import matplotlib.pyplot as plt

ops = eikonal(0.5, 0.5, 0.5, 11, 13, 6)
cu_opt = ops.solve()


def distToHyp(X,Y):
  xmy = X-Y
  hptX = 1/((Y/X)+1)
  hptY = 1 - 1/((Y/X)+1)
  theta = np.arccos(xmy / np.sqrt(2*(X**2 + Y**2)))
  dh = np.sin(theta) * np.sqrt((X-hptX)**2 + (Y-hptY)**2)
  return dh
  

def distToTri(X,Y):
  Z = distToHyp(X,Y)
  mins = np.zeros_like(X)
  N = np.size(X,0)
  for j in range(N):
    mins[j] = np.min([X[j],Y[j],Z[j]])
  return mins


np.savetxt("cu_opt_new.txt", cu_opt)
#print(cu_opt)
#
#v00 = ops.evalPoly(0.25,0.25)
#print(v00.dot(cu_opt))

Zfine = np.loadtxt("triquadleg_n22_m35_N253.txt")
Nfine = 253; 
Xfine = Zfine[0:Nfine]; 
Yfine = Zfine[Nfine:2*Nfine] 
jP = jPoly(Nfine, ops.n, ops.a, ops.b, ops.c, 4)
Vfine = computeV(jP, ops.N, Xfine, Yfine)
tmins = Vfine.dot(cu_opt); #tmins = tmins / np.max(tmins)
dmins = distToTri(Xfine, Yfine); dmins = dmins / np.max(dmins)


fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
#ax.plot_trisurf(Xfine, Yfine, tmins)
ax.plot_trisurf(Xfine, Yfine, tmins)
#ax.scatter(Xfine, Yfine, dmins)
plt.show()


