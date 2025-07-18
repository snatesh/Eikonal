import numpy as np


speed_of_sound = 34300.0 # cm/s

class SimulatedUltrasonicSensor:
  def __init__(self, x, y, theta, obstacles, room_dims, max_range=100, step=0.5):
    """
    Initialize a simulated ultrasonic sensor.
    x, y: position of the sensor in world coordinates (cm)
    theta: orientation of the sensor in radians
    obstacles: list of world-space obstacles (used for raycasting)
    """
    self.x = x
    self.y = y
    self.theta = theta
    self.obstacles = obstacles
    self.room_dims = room_dims
    self.max_range = max_range
    self.step = step;
  
  def trigger_and_echo_time(self):
    """
    Simulate sending a pulse and waiting for echo.
    Returns the round-trip time in seconds.
    """
    return simulate_sensor(self.x, self.y, self.theta, \
                           self.obstacles, self.max_range, self.step)
  
  def get_distance_cm(self):
    """
    Convert round-trip time-of-flight to a one-way distance in cm.
    """
    tof = self.trigger_and_echo_time()
    return (tof * speed_of_sound) / 2

def transform_sensor_to_world(robot_pose, sensor_rel):
  """
  Transform sensor pose from robot frame to world frame.

  Parameters:
    robot_pose: (x, y, theta) of robot in world frame (theta in radians)
    sensor_rel: (dx, dy, dtheta) of sensor in robot frame (dtheta in degrees)

  Returns:
    (x_s, y_s, theta_s): sensor pose in world frame
  """
  x, y, theta = robot_pose
  dx, dy, dtheta = sensor_rel
  dtheta = dtheta

  cos_t = np.cos(theta)
  sin_t = np.sin(theta)

  x_s = x + cos_t * dx - sin_t * dy
  y_s = y + sin_t * dx + cos_t * dy
  theta_s = theta + dtheta

  return x_s, y_s, theta_s

def simulate_sensor(x, y, theta, obstacles, max_range, step, cone_deg=35, n_rays=7):
  """
  Simulate ultrasonic sensor reading with a conical beam pattern.
  Casts multiple rays within a cone centered at `theta`.

  Parameters:
    x, y        : sensor world coordinates (cm)
    theta       : sensor central direction (radians)
    obstacles   : list of axis-aligned boxes ((xmin, ymin), (xmax, ymax))
    max_range   : maximum detection distance (cm)
    step        : distance increment for raymarching (cm)
    cone_deg    : total cone angle in degrees
    n_rays      : number of rays to cast within the cone

  Returns:
    shortest round-trip time (seconds) among all rays
  """
  cone_rad = np.radians(cone_deg)
  half_angle = cone_rad / 2
  angles = np.linspace(theta - half_angle, theta + half_angle, n_rays)

  min_d = max_range
  for angle in angles:
    for d in np.arange(0, max_range, step):
      px = x + d * np.cos(angle)
      py = y + d * np.sin(angle)

      for (xmin, ymin), (xmax, ymax) in obstacles:
        if xmin <= px <= xmax and ymin <= py <= ymax:
          min_d = min(min_d, d)
          break  # hit found for this ray
      else:
        continue  # no hit, keep stepping this ray
      break  # break outer loop on first hit for this ray

  return (2 * min_d) / speed_of_sound  # round-trip time in seconds


def setup_sensors(robot_pose, sensor_offsets, obstacles, room_dims, max_range, step):
  """
  Given the robot's pose and sensor mounting offsets,
  create and return a list of SimulatedUltrasonicSensor objects.
  """
  sensors = []
  for sensor_rel in sensor_offsets:
    x_s, y_s, theta_s = transform_sensor_to_world(robot_pose, sensor_rel)
    sensor = SimulatedUltrasonicSensor(x_s, y_s, theta_s, obstacles, room_dims, max_range, step)
    sensors.append(sensor)
  return sensors


