import nlopt
from numpy import *
from eikonal import *
import matplotlib.pyplot as plt

ops = eikonal(0.5, 0.5, 0.5, 11, 13)

def f(cu, grad=None):
  return ops.F(cu)

def cl(result, cu, grad=None):
  result[:] = ops.vl.dot(cu)

def cb(result, cu, grad=None):
  result[:] = ops.vb.dot(cu)

def ch(result, cu, grad=None):
  result[:] = ops.vh.dot(cu)


opt = nlopt.opt(nlopt.LN_SBPLX, ops.N)
print(opt.get_algorithm_name())
opt.set_min_objective(f)
tol = np.ones((ops.N,)) * 1e-14
opt.set_lower_bounds(-np.inf)
opt.set_upper_bounds(np.inf)
#opt.add_equality_mconstraint(cl, tol)
#opt.add_equality_mconstraint(cb, tol)
#opt.add_equality_mconstraint(ch, tol)

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
Nfine = 253; Xfine = Zfine[0:Nfine]; Yfine = Zfine[Nfine:2*Nfine] 
jP = jPoly(Nfine, ops.n, ops.a, ops.b, ops.c, 6)
Vfine = computeV(jP, ops.N, Xfine,Yfine)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.plot(Xfine,Yfine,'k.')
ax.scatter(Xfine, Yfine, np.abs(Vfine.dot(cu_opt)))
#ax.scatter(ops.X, ops.Y, ops.V.dot(cu_opt))
plt.show()


