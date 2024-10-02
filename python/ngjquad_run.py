from ngjQuad import *

n = 23
m = 35
a = 0.5
b = 0.5
c = 0.5
tol = 1e-14
tolc = 1e-14
alph = 0
use_newton = False
use_wolfe = False
nthreads = 6

Z = ngjQuad ( n, m, a, b, c, tol, tolc, 
               alph, use_newton, use_wolfe,
               nthreads )

np.savetxt("high_order_quad", Z)

