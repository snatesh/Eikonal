import nlopt
from numpy import *
from eikonal import *
import matplotlib.pyplot as plt

ops = eikonal(0.5, 0.5, 0.5, 11, 13, 6)
cu_opt = np.loadtxt("cu_opt_new.txt")

def f(cu, grad=None):
  return ops.F(cu)

def cl(result, cu, grad=None):
  result[:] = ops.vl.dot(cu)

def cb(result, cu, grad=None):
  result[:] = ops.vb.dot(cu)

def ch(result, cu, grad=None):
  result[:] = ops.vh.dot(cu)

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

def tracePathLocal(ops, cu, x, y):
  X = ops.Xcirc + x
  Y = ops.Ycirc + y
  V = computeV(ops.polyCirc, ops.N, X, Y) 
  t = V.dot(cu)
  minInd = np.argmin(t)
  return np.array([X[minInd], Y[minInd], t[minInd]])

def tracePathGlobal(ops, cu, x, y):
  if x == y:
    x = x + 1e-3
    y = y - 1e-3 
  path = []
  path.append(np.array([x,y]))
  minX = tracePathLocal(ops, cu_opt, x, y)
  path.append(np.array([minX[0], minX[1]]))
  while minX[2] > 1e-5:
    minX = tracePathLocal(ops, cu_opt, minX[0], minX[1])
    path.append(np.array([minX[0], minX[1]]))
    if  minX[0] < 0 or minX[1] < 0 or\
        minX[0] > 1 or minX[1] > 1 or\
        minX[1] > 1 - minX[0]:
      break
  return np.ndarray((len(path), 2), buffer = np.array(path))



Zfine = np.loadtxt("triquadleg_n22_m35_N253.txt")
Nfine = 253; 
Xfine = Zfine[0:Nfine]; 
Yfine = Zfine[Nfine:2*Nfine] 

xx, yy = np.meshgrid(np.linspace(0,1,20), np.linspace(0,1,20))

X = xx[yy[:]<1-xx[:]]
Y = yy[yy[:]<1-xx[:]]



for j in range(np.size(X)):
  path = tracePathGlobal(ops, cu_opt, X[j], Y[j])
  plt.plot(path[:,0], path[:,1],'.-')
  

plt.plot([0,1],[0,0],'k-')
plt.plot([0,0],[0,1],'k-')
plt.plot([0,1],[1,0],'k-')
#plt.plot(Xfine, Yfine, 'o')
plt.plot(X, Y,'s')
plt.show()


#v00 = ops.evalPoly(0.25,0.25)
#print(v00.dot(cu_opt))
#
##x = 0.5; y = 0.5; 
#thetas = np.linspace(0,2*np.pi, 10)
#path = tracePathGlobal(cu_opt, 0.25, 0.25, 0.001, thetas)
#print(path)
#

#jP = jPoly(Nfine, ops.n, ops.a, ops.b, ops.c, 4)
#Vfine = computeV(jP, ops.N, Xfine, Yfine)
#tmins = Vfine.dot(cu_opt); tmins = tmins / np.max(tmins)
#dmins = distToTri(Xfine, Yfine); dmins = dmins / np.max(dmins)
#
#
#fig = plt.figure(1)
#ax = fig.add_subplot(111, projection='3d')
##ax.plot_trisurf(Xfine, Yfine, tmins)
#ax.scatter(Xfine, Yfine, tmins)
#ax.scatter(Xfine, Yfine, dmins)
#plt.show()
#
#
