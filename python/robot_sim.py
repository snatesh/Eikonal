import math
import numpy as np
from SimulatedUltrasonicSensor import *
from RobotSimConfig import *
from generate_omegaR import *
import matplotlib.pyplot as plt

Rq = np.loadtxt("../bin/xtri_N496_n30_M1378_m51.txt");
Sq = np.loadtxt("../bin/ytri_N496_n30_M1378_m51.txt");
Wq = np.loadtxt("../bin/wtri_N496_n30_M1378_m51.txt"); 

print(Rq.shape)

# max range of sensor (cm)
max_range = 150
# step size for ray casting (cm)
ray_step = 0.2
# Define room and obstacles
room = (400, 400)  # in cm
#obstacles = [((30, 30), (40, 60)), ((60, 20), (70, 50))]  # list of rectangles (xmin, ymin), (xmax, ymax)
obstacles = generate_random_obstacles(10, room)

# Robot config
robot_pose = (50, 50, 0)  # x, y, theta (in rad)
half_size = 5.5           # half robot size (cm)
sensor_offsets = [        # rel sensor positions (x,y,angle)
  (5.5,  3.0,  1.1781),   # front-left outer
  (5.5,  1.0,  0.3927),   # front-left inner
  (5.5, -1.0, -0.3927),   # front-right inner
  (5.5, -3.0, -1.1781)    # front-right outer
]                         # robot is square ~11cm^2 

# params 
cone_angle_deg = 165 # cone spread in degrees
n_sweep = 12  # e.g. sweep across cone in 7 angular steps
dtheta = np.radians(cone_angle_deg / (n_sweep - 1))  # per-rotation step
# angles to simulate robot in-place rotation for sonic sweep
sweep_angles = np.linspace(-cone_angle_deg / 2, cone_angle_deg / 2, n_sweep)
sweep_angles = np.radians(sweep_angles)  # convert to radians

n_steps = 1000
max_step = 2.0 # cm
max_turn = np.pi/12.0 # 15 deg

robot_path = [robot_pose]
fig, ax = plt.subplots(figsize=(8, 8))
for _ in range(n_steps):
  
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
  sector, tris = generate_sector_tris(robot_pose, cone_angle_deg, max_range, n_points=10)
  # quad points mapped to each triangle in sector
  XX, YY = map_quad_grid(Rq, Sq, tris)
  
  
  # Draw sector polygon
  ax.plot(sector[:, 0], sector[:, 1], 'r-', label='Sensor Sector Boundary')
  ax.plot([sector[-1,0],sector[0,0]], [sector[-1,1],sector[0,1]],'r-')
  ax.fill(sector[:, 0], sector[:, 1], color='red', alpha=0.1)
  
  
  
  # Plot robot, sensors cones and splatter
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
  
  
  hit_points = np.array(hit_points)
   
  if hit_points.size != 0: 
    ax.plot(hit_points[:,0], hit_points[:,1], 'cx', label='Hit Points')
  
  rho_vals = rho_splatter(XX, YY, hit_points, sensor_points)
  speed = 1 - rho_vals + 1e-6
  
  
  c = ax.tricontourf(XX.ravel(), YY.ravel(), speed.ravel(), levels=30, cmap='hot')
  colorbar = fig.colorbar(c, ax=ax, label='f(x, y)') 
  
  # Plot obstacles
  for (xmin, ymin), (xmax, ymax) in obstacles:
    rect = plt.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, color='gray', alpha=0.4)
    ax.add_patch(rect)
  
  plt.legend()
  plt.pause(0.1)
  colorbar.remove()
  ax.cla()
  x, y, theta_new = update_pose(room, robot_pose, obstacles, max_step, max_turn, half_size)
  robot_pose = (x, y, theta_new)
  robot_path.append(robot_pose)




