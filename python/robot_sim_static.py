import math
import numpy as np
from SimulatedUltrasonicSensor import *
from RobotSimConfig import *
from generate_omegaR import *
import matplotlib.pyplot as plt

Rq = np.loadtxt("../bin/xtri_N496_n30_M1378_m51.txt");
Sq = np.loadtxt("../bin/ytri_N496_n30_M1378_m51.txt");
Wq = np.loadtxt("../bin/wtri_N496_n30_M1378_m51.txt"); 

# max range of sensor (cm)
max_range = 150
# step size for ray casting (cm)
ray_step = 0.2
# Define room and obstacles
room = (400, 400)  # in cm
#obstacles = [((100, 0), (150, 150))]  # list of rectangles (xmin, ymin), (xmax, ymax)
obstacles = generate_random_obstacles(20, room)

# Robot config
#robot_pose = (50, 50, 0)  # x, y, theta (in rad)
robot_pose = (179, 192, 0)#np.pi/4)  # x, y, theta (in rad)
half_size = 5.5           # half robot size (cm)
# sensor config
sensor_offsets = [        # rel sensor positions (x,y,angle)
  (5.5,  3.0,  1.1781),   # front-left outer
  (5.5,  1.0,  0.3927),   # front-left inner
  (5.5, -1.0, -0.3927),   # front-right inner
  (5.5, -3.0, -1.1781)    # front-right outer
]                         # robot is square ~11cm^2 
cone_angle_deg = 165
n_sweep = 12  # e.g. sweep across cone in 7 angular steps
# angles to simulate robot in-place rotation for sonic sweep
theta_start = robot_pose[2] + np.pi/2.0
theta_end   = robot_pose[2] - np.pi/2.0
sweep_angles = np.linspace(theta_start, theta_end, n_sweep)

# Plot setup
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.set_xlim(0, room[0])
ax.set_ylim(0, room[0])
ax.set_title("AlphaBot2 Sensor Coverage")
ax.set_xlabel("x (cm)")
ax.set_ylabel("y (cm)")
ax.grid(True)

# Plot robot body (11×11 cm square)
robot_x, robot_y, robot_theta = robot_pose
robot_box = np.array([
  [-half_size, -half_size],
  [ half_size, -half_size],
  [ half_size,  half_size],
  [-half_size,  half_size],
  [-half_size, -half_size]
])
R = np.array([
  [np.cos(robot_theta), -np.sin(robot_theta)],
  [np.sin(robot_theta),  np.cos(robot_theta)]
])

robot_box_world = (R @ robot_box.T).T + np.array([robot_x, robot_y])
ax.plot(robot_box_world[:, 0], robot_box_world[:, 1], 'b-', label="Robot")

# Generate sector (Omega_R) a conical sector spanning from robot center
# and also compute triangulation
#sector, tris = generate_sector_tris(robot_pose, cone_angle_deg, max_range, n_points=10)
OmegaR, tris, grid_world, A, b = generate_sweep_rectangle(robot_pose, cone_angle_deg, max_range, 5)

# quad points mapped to each triangle in sector
XX, YY = map_quad_grid(Rq, Sq, tris)

# Draw rectangular OmegaR
ax.plot(OmegaR[:, 0], OmegaR[:, 1], 'r-', label='Sensor Sector Boundary')
ax.plot([OmegaR[-1,0],OmegaR[0,0]], [OmegaR[-1,1],OmegaR[0,1]],'r-')
ax.fill(OmegaR[:, 0], OmegaR[:, 1], color='red', alpha=0.1)
## Plot robot, sensors cones and splatter
hit_points = []
sensor_points = []

for dth in sweep_angles:
  # temporarily adjust robot heading
  theta_sweep = robot_pose[2] + dth
  robot_pose_sweep = (robot_pose[0], robot_pose[1], theta_sweep) 
  # simulate sensor readings at this pose
  for i, sensor_rel in enumerate(sensor_offsets):
    x_s, y_s, theta_s = transform_sensor_to_world(robot_pose_sweep, sensor_rel)

    d_cm, hit_pt = simulate_single_beam_sensor(x_s, y_s, theta_s, obstacles, max_range, ray_step)

    if hit_pt is not None:
      hit_points.append(hit_pt)
      sensor_points.append((x_s, y_s))


OmegaR_c = (A@(OmegaR-OmegaR[0]).T + b[:,None]).T
fig1, ax1 = plt.subplots(figsize=(8, 8))
ax1.plot(OmegaR_c[:, 0], OmegaR_c[:, 1], 'r-', label='Sensor Sector Boundary')
ax1.plot([OmegaR_c[-1,0],OmegaR_c[0,0]], [OmegaR_c[-1,1],OmegaR_c[0,1]],'r-')


hit_points = np.array(hit_points)
print(hit_points.shape)
sensor_points = np.array(sensor_points)

ax.plot(hit_points[:,0], hit_points[:,1], 'cx', label='Hit Points')

hit_points_can = (A @ (hit_points - OmegaR[0]).T + b[:,None]).T
sensor_points_can = (A @ (sensor_points - OmegaR[0]).T + b[:,None]).T


rho_vals = rho_splatter(XX, YY, hit_points_can, sensor_points_can)
speed = 1 - rho_vals + 1e-6

ax1.plot(hit_points_can[:,0], hit_points_can[:,1], 'cx', label='Hit Points')
c = ax1.tricontourf(XX.ravel(), YY.ravel(), speed.ravel(), levels=30, cmap='hot')
colorbar = fig1.colorbar(c, ax=ax1, label='f(x, y)') 
ax1.plot(tris.points[:, 0], tris.points[:, 1], 'o', label=None) # Plot the points
ax1.triplot(tris.points[:, 0], tris.points[:, 1], tris.simplices, label='Delaunay Triangulation') # Plot the triangles


# Plot obstacles
for (xmin, ymin), (xmax, ymax) in obstacles:
  rect = plt.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, color='gray', alpha=0.4)
  ax.add_patch(rect)

np.savetxt("tri_points.txt", tris.points, fmt="%18e")
np.savetxt("tri_faces.txt", tris.simplices, fmt="%d")
np.savetxt("speed_field.txt", speed, fmt="%.18e")

plt.legend()
plt.show()





