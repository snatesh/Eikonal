from jPoly import *
import matplotlib.pyplot as plt

Qrule = np.loadtxt("triquadLeg_10_17.txt"); 
N = int(len(Qrule) / 3)
X = Qrule[0:N]
Y = Qrule[N:2*N]
W = Qrule[2*N:3*N]

g = lambda x, y : x**2 + y**2
cg, V = computeCoeffs(g, 17, 0.5, 0.5, 0.5, X, Y, W, 1)
ga = V.dot(cg)
print(ga.shape)
print(cg)
print(g(X,Y).shape)
print(np.linalg.norm(ga - g(X,Y)))

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
#ax.plot_trisurf(X, Y, g(X,Y))
ax.plot_trisurf(X, Y, ga)
plt.show()


