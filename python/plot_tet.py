import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt('../bin/xtet.txt') * 6.0;
y = np.loadtxt('../bin/ytet.txt') * 6.0;
z = np.loadtxt('../bin/ztet.txt') * 6.0;


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

