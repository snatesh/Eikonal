import numpy as np
import matplotlib.pyplot as plt

dir = "../testing/testdata/"

x = np.loadtxt(dir + 'xtet_N35_n4_M35_m5.txt');
y = np.loadtxt(dir + 'ytet_N35_n4_M35_m5.txt');
z = np.loadtxt(dir + 'ztet_N35_n4_M35_m5.txt');



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

