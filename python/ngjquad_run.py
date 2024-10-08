from ngjQuad import *

n = 22
m = 26
a = 0.5
b = 0.5
c = 0.5
tol = 1e-14
tolc = 1e-14
alph = 0
use_newton = False
use_wolfe = False
nthreads = 6
N = int(0.5*n*(n+1))
fname = 'triquadleg_n' + str(n) + '_m' \
          + str(m) + '_N' + str(N) + '.txt'
print(fname)
Z = ngjQuad ( n, m, a, b, c, tol, tolc, 
               alph, use_newton, use_wolfe,
               nthreads )


np.savetxt(fname, Z)
