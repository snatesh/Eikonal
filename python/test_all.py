from pMat import *
from sFactors import *
from jPoly import *
import matplotlib.pyplot as plt

def distToHyp(X,Y):
  Np = np.size(X,0)
  dh = np.zeros_like(X)
  for j in range(Np):
    x = np.array([X[j],Y[j]])
    z = np.array([(X[j]-Y[j] + 1) /2,
                  (Y[j]-X[j] + 1) / 2])
    dh[j] = np.linalg.norm(x - z)
  return dh


def f(X,Y):
  Z = distToHyp(X,Y)
  mins = np.zeros_like(X)
  N = np.size(X,0)
  for j in range(N):
    mins[j] = np.min([X[j],Y[j],Z[j]])
  return mins

#def f(X,Y):
#  return np.sin(X**2 + Y**2)

# jacobi poly params
a = 0.5; b = 0.5; c = 0.5; n = 20
Z = np.loadtxt("triquadLeg_16_27.txt")
N = int(len(Z) / 3.0)
X = Z[0:N]
Y = Z[N:2*N]
W = Z[2*N:3*N]
poly = jPoly(N, n, a, b, c, 1)
c, V = computeCoeffs(poly, f, X, Y, W)
print(c)
print(W.dot(f(X,Y))/2)
print(W.dot(V.dot(c))/2)
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, Y, V.dot(c), alpha = 0.7)
ax.plot_trisurf(X, Y, f(X,Y), alpha = 0.4)
plt.show()

