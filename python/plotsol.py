import matplotlib.pyplot as plt
import numpy as np

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


Z = np.loadtxt("triquadLeg_16_27.txt");
U = np.loadtxt("Usol.txt");
U = U / np.max(U)

N = int(len(Z) / 3.0)

X = Z[0:N]
Y = Z[N:2*N]

Uref = distToTri(X, Y)
Uref = Uref / np.max(Uref)

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, Y, U, alpha = 0.8)
ax.plot_trisurf(X, Y, Uref, alpha = 0.7)
plt.show()

