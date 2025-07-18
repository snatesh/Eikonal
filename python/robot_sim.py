import math
import numpy as np
from SimulatedUltrasonicSensor import *
from RobotSimConfig import *
import matplotlib.pyplot as plt

# max range of sensor (cm)
max_range = 100
# step size for ray casting (cm)
ray_step = 0.2
# Define room and obstacles
room = (100, 100)  # in cm
#obstacles = [((30, 30), (40, 60)), ((60, 20), (70, 50))]  # list of rectangles (xmin, ymin), (xmax, ymax)
obstacles = generate_random_obstacles(4, room, min_size=5, max_size=15, seed=None)

# Robot config
robot_pose = (50, 50, 0)  # x, y, theta (in rad)
half_size = 5.5           # half robot size (cm)
sensor_offsets = [        # rel sensor positions (x,y,angle)
  (5.5,  3.0,  1.1781),   # front-left outer
  (5.5,  1.0,  0.3927),   # front-left inner
  (5.5, -1.0, -0.3927),   # front-right inner
  (5.5, -3.0, -1.1781)    # front-right outer
]                         # robot is square ~11cm^2 

n_steps = 1000
max_step = 2.0 # cm
max_turn = np.pi/12.0 # 15 deg

robot_path = [robot_pose]
fig, ax = plt.subplots(figsize=(8, 8))
for _ in range(n_steps):
  
  for i, sensor_rel in enumerate(sensor_offsets):
    x_s, y_s, theta_s = transform_sensor_to_world(robot_pose, sensor_rel)
    print(f"Sensor {i}: x = {x_s:.2f}, y = {y_s:.2f}, theta = {theta_s:.2f}")
  
  
  simulated_sensors = setup_sensors(robot_pose, sensor_offsets, obstacles, room, max_range, ray_step)
  distances = [sensor.get_distance_cm() for sensor in simulated_sensors]
  print("Sensor readings (cm):", distances)
  
  # Plot setup

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
  
  # Plot obstacles
  for (xmin, ymin), (xmax, ymax) in obstacles:
    rect = plt.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, color='gray', alpha=0.4)
    ax.add_patch(rect)
  
  # Plot sensors and cones
  cone_deg = 35
  cone_length = 25  # cm
  for i, (sensor_rel, reading_cm) in enumerate(zip(sensor_offsets, distances)):
    x_s, y_s, theta_s = transform_sensor_to_world(robot_pose, sensor_rel)
    (line,) = ax.plot(x_s, y_s, 'o'); clr = line.get_color()
    ax.text(x_s-2, y_s-1, f"S{i}", color='k', ha='center', fontsize=10)
  
    half_angle = np.radians(cone_deg / 2)
    for angle in [theta_s - half_angle, theta_s, theta_s + half_angle]:
      x_end = x_s + cone_length * np.cos(angle)
      y_end = y_s + cone_length * np.sin(angle)
      slabel = f"Sensor {i}: {reading_cm:.1f} cm" if angle == theta_s else None
      ax.plot([x_s, x_end], [y_s, y_end], '--', linewidth=1, color = clr, label=slabel)
  
  plt.legend()
  plt.pause(1)
  ax.cla()
  x, y, theta_new = update_pose(room, robot_pose, obstacles, max_step, max_turn, half_size)
  robot_pose = (x, y, theta_new)
  robot_path.append(robot_pose)




