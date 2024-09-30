import nlopt
from numpy import *
from eikonal import *
import matplotlib.pyplot as plt

ops = eikonal(0.5, 0.5, 0.5, 13, 14, 4)

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

opt = nlopt.opt(nlopt.LN_COBYLA, ops.N)
print(opt.get_algorithm_name())
opt.set_min_objective(f)
tol = np.ones((ops.N,)) * 1e-14
opt.set_lower_bounds(-np.inf)
opt.set_upper_bounds(np.inf)
opt.add_equality_mconstraint(cl, tol)
opt.add_equality_mconstraint(cb, tol)
opt.add_equality_mconstraint(ch, tol)

opt.set_ftol_rel(1e-14)
opt.set_xtol_abs(1e-14)
cu = np.zeros((ops.N,))
cu[0] = 1
#np.random.randn(ops.N)
cu_opt = opt.optimize(cu)
print(cu_opt)

v00 = ops.evalPoly(0.1,0.1)
print(v00.dot(cu_opt))

Zfine = np.loadtxt("triquadleg_n22_m35_N253.txt")
Nfine = 253; 
Xfine = Zfine[0:Nfine]; 
Yfine = Zfine[Nfine:2*Nfine] 
jP = jPoly(Nfine, ops.n, ops.a, ops.b, ops.c, 4)
Vfine = computeV(jP, ops.N, Xfine, Yfine)
tmins = Vfine.dot(cu_opt); tmins = tmins / np.max(tmins)
dmins = distToTri(Xfine, Yfine); dmins = dmins / np.max(dmins)


fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
#ax.plot_trisurf(Xfine, Yfine, tmins)
ax.scatter(Xfine, Yfine, tmins)
ax.scatter(Xfine, Yfine, dmins)
plt.show()


