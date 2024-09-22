import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt('../testing/testdata/xtet_56.txt') * 6.0;
y = np.loadtxt('../testing/testdata/ytet_56.txt') * 6.0;
z = np.loadtxt('../testing/testdata/ztet_56.txt') * 6.0;


fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(projection='3d')
ax.scatter(x, y, z, c = 'r', s = 50)
tet_e1 = np.linspace(0,1);
ax.plot(tet_e1, np.zeros_like(tet_e1), np.zeros_like(tet_e1), 'k-', linewidth=3)
ax.plot(np.zeros_like(tet_e1), tet_e1, np.zeros_like(tet_e1), 'k-', linewidth=3)
ax.plot(np.zeros_like(tet_e1), np.zeros_like(tet_e1), tet_e1, 'k-', linewidth=3)
ax.plot(tet_e1, (1-tet_e1), np.zeros_like(tet_e1), 'k-', linewidth=3)
ax.plot(tet_e1, np.zeros_like(tet_e1), (1-tet_e1), 'k-', linewidth=3)
ax.plot(np.zeros_like(tet_e1), tet_e1, (1-tet_e1), 'k-', linewidth=3)
plt.show()

