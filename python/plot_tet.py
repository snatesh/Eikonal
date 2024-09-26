import numpy as np
import matplotlib.pyplot as plt

dir = "../testing/testdata/"


#x = np.loadtxt(dir + 'xtet_N56_n5_M120_m7.txt');
#y = np.loadtxt(dir + 'ytet_N56_n5_M120_m7.txt');
#z = np.loadtxt(dir + 'ztet_N56_n5_M120_m7.txt');

x = np.loadtxt("../bin/xtet_new.txt");
y = np.loadtxt("../bin/ytet_new.txt");
z = np.loadtxt("../bin/ztet_new.txt");
w = np.loadtxt("../bin/Wtet.txt");

print(np.sum(w));

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

