import nlopt
from numpy import *
from eikonal import *
import matplotlib.pyplot as plt

nthreads = 6

def fzero(X,Y):
  return np.zeros_like(X)

ops = eikonal(fzero, 0.5, 0.5, 0.5, 14, 16, nthreads)
cu_opt = np.loadtxt("cu_opt_new1.txt")



def tracePathLocal(ops, cu, x, y):
  X = ops.Xcirc + x
  Y = ops.Ycirc + y
  V = computeV(ops.polyCirc, ops.N, X, Y) 
  t = V.dot(cu)
  minInd = np.argmin((t))
  return np.array([X[minInd], Y[minInd], (t[minInd])])

def tracePathGlobal(ops, cu, x, y):
  path = []
  path.append(np.array([x,y]))
  minX = tracePathLocal(ops, cu_opt, x, y)
  path.append(np.array([minX[0], minX[1]]))
  buf = 1e-3
  while minX[2] > 1e-5:
    minX = tracePathLocal(ops, cu_opt, minX[0], minX[1])
    path.append(np.array([minX[0], minX[1]]))
    if  minX[0] < buf or minX[1] < buf or\
        minX[0] > 1-buf or minX[1] > 1-buf or\
        minX[1] > (1 - minX[0])-buf:
      break
  return np.ndarray((len(path), 2), buffer = np.array(path))

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

#Zfine = np.loadtxt("triquadleg_n22_m35_N253.txt")
#Nfine = 253; 
#Xfine = Zfine[0:Nfine]; 
#Yfine = Zfine[Nfine:2*Nfine] 
jP = jPoly(ops.N, ops.n, ops.a+1, ops.b+1, ops.c+1, nthreads, True)
V = computeV(jP, ops.N, ops.X, ops.Y)
ops.W = ops.W / 2.0
print(np.sum(ops.W))
print(np.linalg.norm(ops.W.dot(V.dot(cu_opt) - distToTri(ops.X,ops.Y)))/np.linalg.norm(ops.W.dot(distToTri(ops.X,ops.Y))))
xx, yy = np.meshgrid(np.linspace(0,1,20), np.linspace(0,1,20))

X = xx[yy[:]<1-xx[:]]
Y = yy[yy[:]<1-xx[:]]
#nfine = np.size(X,0)
#jP = jPoly(nfine, ops.n, ops.a, ops.b, ops.c, 6)
#Vfine = computeV(jP, ops.N, X, Y)
#tmins = Vfine.dot(cu_opt); tmins = tmins / np.max(tmins)
#dmins = distToTri(X, Y); dmins = dmins / np.max(dmins)

for j in range(np.size(X)):
  path = tracePathGlobal(ops, cu_opt, X[j], Y[j])
  plt.plot(path[:,0], path[:,1],'.-')

plt.plot([0,1],[0,0],'k-')
plt.plot([0,0],[0,1],'k-')
plt.plot([0,1],[1,0],'k-')
plt.plot(X, Y,'s')
plt.show()

##
#fig = plt.figure(1)
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_trisurf(X, Y, tmins)
#ax.plot_trisurf(X, Y, dmins)
#plt.show()
##
##
