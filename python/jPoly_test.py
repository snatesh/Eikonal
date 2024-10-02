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
  

def distToTri(X,Y):
  Z = distToHyp(X,Y)
  mins = np.zeros_like(X)
  N = np.size(X,0)
  for j in range(N):
    mins[j] = np.min([X[j],Y[j],Z[j]])
  return mins

def finv(X, Y):
  return X*Y*(1-X-Y)
  #centX = 1./3.
  #centY = 1./3.
  #return np.exp(-( (X-centX)**2 + (Y-centY)**2 )) * X * Y
  #F = np.zeros_like(dist)
  #for j  in range(np.size(dist,0)):
  #  if dist[j] < 1e-2:
  #    F[j] = 20
  #  else:
  #    F[j] = 1
  #return 1.0 / F  



Qrule = np.loadtxt("triquadLeg_10_17.txt"); 
N = int(len(Qrule) / 3)
X = Qrule[0:N]
Y = Qrule[N:2*N]
W = Qrule[2*N:3*N]

cF, V = computeCoeffs(finv, 17, 0.5, 0.5, 0.5, X, Y, W, 1)
Fa = V.dot(cF)
print(np.linalg.norm(Fa - finv(X,Y)))
#
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
#ax.plot_trisurf(X, Y, finv(X,Y))
ax.plot_trisurf(X, Y, Fa)
plt.show()


