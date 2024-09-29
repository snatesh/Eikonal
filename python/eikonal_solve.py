import nlopt
from numpy import *
from eikonal import *


ops = eikonal(0.5,0.5,0.5)

def f(cu, grad=None):
  return ops.F(cu)


opt = nlopt.opt(nlopt.LN_SBPLX, ops.N)
print(opt.get_algorithm_name())
print(opt.get_dimension())
opt.set_min_objective(f)
opt.set_lower_bounds(-np.inf)
opt.set_upper_bounds(np.inf)

opt.set_ftol_rel(1e-14)
opt.set_xtol_abs(1e-14)


cu = np.random.randn(ops.N)
cu_opt = opt.optimize(cu)
print(cu)
print(cu_opt)
