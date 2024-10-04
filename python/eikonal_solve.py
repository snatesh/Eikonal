from numpy import *
from eikonal import *
import matplotlib.pyplot as plt

def fzero(X,Y):
  return np.zeros_like(X)

def finv(X, Y):
  return np.asfortranarray(np.ones_like(X))
  #return np.asfortranarray(X*Y*(1-X-Y))
  #centX = 1./3.
  #centY = 1./3.
  #return 1./(np.exp(-( (X-centX)**2 + (Y-centY)**2 )) * X * Y)

#opsH = eikonal(finv, 0.5, 0.5, 0.5, 14, 16, 6)
#cu_optH = opsH.solveH()
#ul = opsH.vl.dot(cu_optH)
#ub = opsH.vb.dot(cu_optH)
#uh = opsH.vh.dot(cu_optH)
#opsP = eikonal(fzero, 0.5, 0.5, 0.5, 14, 16, 6)
ops = eikonal(finv, 0.5, 0.5, 0.5, 14, 16, 6)
ul = np.zeros((ops.N,))#opsP.vl.dot(cu_optH)
ub = np.zeros((ops.N,))#opsP.vb.dot(cu_optH)
uh = np.zeros((ops.N,))#opsP.vh.dot(cu_optH)

cu_opt = ops.solveP(ul, ub, uh)

#cu_opt = cu_optH + cu_optP

def distToHyp(X,Y):
  xmy = X-Y
  hptX = 1/((Y/X)+1)
  hptY = 1 - 1/((Y/X)+1)
  theta = np.arccos(xmy / np.sqrt(2*(X**2 + Y**2)))
  dh = np.sin(theta) * np.sqrt((X-hptX)**2 + (Y-hptY)**2)
  return dh
  


def u(x1, x2):
  dists = np.zeros_like(x1)
  for j in range(np.size(x1,0)):
    x = np.array([x1[j],x2[j]])
    z = np.array([(x1[j]-x2[j]+1.0)/2.0, (x2[j]-x1[j]+1.0)/2.0])
    d = np.linalg.norm(x-z)
    dists[j] = np.min([x1[j],x2[j],d])
  return dists


np.savetxt("cu_opt_new1.txt", cu_opt)

Zfine = np.loadtxt("triquadleg_n22_m34_N276.txt")
Nfine = 276; 
Xfine = Zfine[0:Nfine]; 
Yfine = Zfine[Nfine:2*Nfine] 
jP = jPoly(Nfine, ops.n, ops.a, ops.b, ops.c, 4)
Vfine = computeV(jP, ops.N, Xfine, Yfine)
tmins = Vfine.dot(cu_opt); tmins = tmins / np.max(tmins)
dmins = u(Xfine, Yfine); dmins = dmins / np.max(dmins)


fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
#ax.plot_trisurf(Xfine, Yfine, tmins)
ax.plot_trisurf(Xfine, Yfine, tmins, alpha = 0.7)
ax.plot_trisurf(Xfine, Yfine, dmins, alpha = 0.7)
#ax.scatter(Xfine, Yfine, dmins)
plt.show()


