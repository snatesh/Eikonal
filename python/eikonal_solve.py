from numpy import *
from eikonal import *
import matplotlib.pyplot as plt

nthreads = 6;

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
ops = eikonal(finv, 0.5, 0.5, 0.5, 14, 16, nthreads)
cu_opt = ops.solveP()

#cu_opt = cu_optH + cu_optP
  
def distToHyp(X,Y):
  Np = np.size(X,0)
  dh = np.zeros_like(X)
  for j in range(Np):
    x = np.array([X[j],Y[j]])
    z = np.array([(X[j]-Y[j] + 1) /2, 
                  (Y[j]-X[j] + 1) / 2])
    dh[j] = np.linalg.norm(x - z)
  return dh
  

def distToTri(X,Y):
  Z = distToHyp(X,Y)
  mins = np.zeros_like(X)
  N = np.size(X,0)
  for j in range(N):
    mins[j] = np.min([X[j],Y[j],Z[j]])
  return mins



np.savetxt("cu_opt_new1.txt", cu_opt)

#Zfine = np.loadtxt("triquadleg_n22_m34_N276.txt")
#Zfine = np.loadtxt("triquadleg_n15_m18_N136.txt")
#Nfine = 276; 
#Nfine = 136; 
#Xfine = Zfine[0:Nfine]; 
#Yfine = Zfine[Nfine:2*Nfine] 
#jP = jPoly(Nfine, ops.n, ops.a, ops.b, ops.c, 6)
#Vfine = computeV(jP, ops.N, Xfine, Yfine)
#tmins = Vfine.dot(cu_opt); tmins = tmins / np.max(tmins)
#dmins = u(Xfine, Yfine); dmins = dmins / np.max(dmins)

xx, yy = np.meshgrid(np.linspace(0,1,40), np.linspace(0,1,40))

X = xx[yy[:]<1-xx[:]]
Y = yy[yy[:]<1-xx[:]]
ngrd = np.size(X,0)
jP = jPoly(ngrd, ops.n, ops.a+1, ops.b+1, ops.c+1, nthreads, False)
V = computeV(jP, ops.N, X, Y)
tmins = V.dot(cu_opt); tmins = tmins / np.max(tmins)
dmins = distToTri(X, Y); dmins = dmins / np.max(dmins)

#fig = plt.figure(1)
#ax = fig.add_subplot(111, projection='3d')
##ax.plot_trisurf(Xfine, Yfine, tmins)
#ax.plot_trisurf(Xfine, Yfine, tmins, alpha = 0.7)
#ax.plot_trisurf(Xfine, Yfine, dmins, alpha = 0.7)
##ax.scatter(Xfine, Yfine, dmins)
#plt.show()
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, Y, tmins, alpha = 0.7)
ax.plot_trisurf(X, Y, dmins, alpha = 0.4)
plt.show()
