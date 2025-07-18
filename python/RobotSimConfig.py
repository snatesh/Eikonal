import numpy as np

def in_obstacle(x, y, obstacles, buffer=5.5):
  """
  Returns True if any part of the robot (with buffered area)
  overlaps any obstacle.
  """
  for (xmin, ymin), (xmax, ymax) in obstacles:
    if (x + buffer > xmin and x - buffer < xmax and
        y + buffer > ymin and y - buffer < ymax):
      return True
  return False

def update_pose(room, pose, obstacles, step_size=2.0, angle_range=np.pi/12, buffer=5.5):
  """
  Updates the robot pose by stepping in a random direction with rotation,
  while staying inside the room and avoiding obstacles.
  """
  x, y, theta = pose
  for _ in range(10):  # Try 10 times
    # Random small rotation
    dtheta = np.random.uniform(-angle_range, angle_range)
    theta_new = theta + dtheta

    # Random step
    step = np.random.uniform(0, step_size)
    dx = step * np.cos(theta_new)
    dy = step * np.sin(theta_new)

    # Clamp position within room bounds with margin
    x_new = np.clip(x + dx, buffer, room[0] - buffer)
    y_new = np.clip(y + dy, buffer, room[1] - buffer)

    if not in_obstacle(x_new, y_new, obstacles, buffer):
      return (x_new, y_new, theta_new)

  # If no valid move after 10 tries, rotate 180 in place
  theta_new = theta + np.pi
  return (x, y, theta_new)

def generate_random_obstacles(n, room_size, min_size=5, max_size=15, seed=None):
  """
  Generate n random rectangular obstacles within the room.
  Each rectangle is axis-aligned, with size between min_size and max_size (in cm).
  """
  if seed is not None:
    np.random.seed(seed)

  width, height = room_size
  obstacles = []

  for _ in range(n):
    w = np.random.uniform(min_size, max_size)
    h = np.random.uniform(min_size, max_size)
    x = np.random.uniform(0, width - w)
    y = np.random.uniform(0, height - h)
    obstacles.append(((x, y), (x + w, y + h)))

  return obstacles
