import numpy as np
from scipy.spatial import Delaunay

def generate_trapezoid(robot_pose, front_width, back_width, height):
  """
  Generates a trapezoid in world coordinates based on robot pose.

  Parameters:
    robot_pose: (x, y, theta) - robot position and heading in radians
    front_width: width of far edge (in cm)
    back_width: width at robot (in cm)
    height: distance forward from robot (in cm)

  Returns:
    np.ndarray of shape (4, 2): world coordinates of trapezoid corners,
    ordered clockwise starting from back-left.
  """
  x, y, theta = robot_pose
  bw2 = back_width / 2
  fw2 = front_width / 2

  # Define corners in robot's local frame
  local_pts = np.array([
    [0, -bw2],         # back-left
    [0,  bw2],         # back-right
    [height, fw2],     # front-right
    [height, -fw2]     # front-left
  ])

  # Rotate and translate into world frame
  R = np.array([
    [np.cos(theta), -np.sin(theta)],
    [np.sin(theta),  np.cos(theta)]
  ])
  world_pts = (R @ local_pts.T).T + np.array([x, y])
  return world_pts


def compute_front_width(cone_angle_deg, height_cm):
  alpha = np.radians(cone_angle_deg)
  return 2 * height_cm * np.tan(alpha / 2)

def generate_sector_tris(robot_pose, cone_angle_deg, max_range_cm, n_points=50):
  """
  Generates a circular sector (wedge) polygon in world coordinates.

  Parameters:
    robot_pose: (x, y, theta) - robot center and heading in radians
    cone_angle_deg: total cone angle (in degrees, e.g. 165)
    max_range_cm: radius of the sector (in cm)
    n_points: number of arc points to define the circular edge

  Returns:
    sector_points: (n_points+2, 2) array of polygon vertices (x, y),
                   starting and ending at robot pose.
  """
  x, y, theta = robot_pose
  half_angle_rad = np.radians(cone_angle_deg / 2)

  # Angles for arc points (relative to heading direction)
  angles = np.linspace(-half_angle_rad, half_angle_rad, n_points)

  # Points along arc in robot-local frame
  arc_points = np.stack([
    max_range_cm * np.cos(angles),
    max_range_cm * np.sin(angles)
  ], axis=1)

  # Rotation matrix for world transform
  R = np.array([
    [np.cos(theta), -np.sin(theta)],
    [np.sin(theta),  np.cos(theta)]
  ])

  # Rotate and translate to world frame
  arc_world = arc_points @ R.T + np.array([x, y])

  # Full polygon = [robot origin] + [arc points] + [robot origin]
  # add robot origin outside if want full sector plot
  # excluded to not have duplicate points for delaunay
  sector_points = np.vstack([
    [x, y],
    arc_world
  ])

  tris = Delaunay(sector_points)

  return sector_points, tris


import numpy as np

def generate_sweep_rectangle(robot_pose, cone_angle_deg, max_range_cm, n):
  """
  Generate a rectangle in front of the robot representing the sensor sweep region.

  Parameters:
    robot_pose        - tuple (x, y, theta) in cm and radians
    cone_angle_deg    - combined cone angle of the sensors (e.g., 165 degrees)
    max_range_cm      - max sensing range in cm

  Returns:
    rect_world: (4,2) ndarray of rectangle corners in world coordinates (CCW)
  """
  x, y, theta = robot_pose
  cone_angle_rad = np.radians(cone_angle_deg)
  #half_width = max_range_cm * np.tan(cone_angle_rad / 2)
  #half_width = max_range_cm / np.tan(np.pi/2-cone_angle_rad)
  half_width = max_range_cm
  # Rectangle in robot frame (origin at robot center)
  # Vertices ordered counter-clockwise
  rect_local = np.array([
    [0, -half_width],      # back-left (robot centerline left)
    [0,  half_width],      # back-right
    [max_range_cm,  half_width],  # front-right
    [max_range_cm, -half_width]   # front-left
  ])

  # Rotate and translate to world frame
  R = np.array([
    [np.cos(theta), -np.sin(theta)],
    [np.sin(theta),  np.cos(theta)]
  ])
  rect_world = (R @ rect_local.T).T + np.array([x, y])

  # map from local to canonical robot frame [0,1]x[-1,1]
  A, b = compute_affine_map_to_canonical(rect_world)

  # Generate grid points in canonical frame
  x_lin = np.linspace(0.0, 1.0, n + 1)
  y_lin = np.linspace(-1.0, 1.0, n + 1)
  XX, YY = np.meshgrid(x_lin, y_lin)
  grid_canonical = np.vstack([XX.ravel(), YY.ravel()]).T  # shape (n*(2n+1), 2)

  # map grid from canonical to local
  grid_local = (A @ (rect_world - rect_world[0]).T + b[:,None]).T 

  # Map grid to world frame
  grid_world = (R @ grid_local.T).T + np.array([x, y])

  # Perform Delaunay triangulation on canonical frame
  tri = Delaunay(grid_canonical)

  return rect_world, tri, grid_world, A, b

def compute_affine_map_to_canonical(rect_world):
  """
  Compute the affine map x_canonical = A @ x_world + b
  such that rect_world → [0,1] × [-1,1]

  Parameters:
    rect_world: (4,2) ndarray of rectangle corners in CCW order
                [bottom-left, bottom-right, top-right, top-left]

  Returns:
    A: (2,2) affine transformation matrix
    b: (2,) translation vector
  """
  # Use origin (bottom-left), x_vec (bottom-right - bottom-left),
  # and y_vec (top-left - bottom-left)
  origin = rect_world[0]
  x_vec = rect_world[3] - origin
  y_vec = rect_world[1] - origin

  # World → canonical map: map x_vec → [1, 0], y_vec → [0, 2]
  # So we solve A @ x_vec = [1, 0], A @ y_vec = [0, 2]
  # Stack: [x_vec, y_vec] is (2,2), so invert and right-multiply
  T_world = np.column_stack((x_vec, y_vec))  # (2,2)
  T_canon = np.array([[1, 0], [0, 2]])       # target vectors

  A = T_canon @ np.linalg.inv(T_world)       # (2,2)
  b = np.array([0,-1])                       # (2,)

  return A, b


def incidence_matrix(Xe):
  """
  Computes the incidence matrix J = [x2 - x1, x3 - x1]
  
  Parameters:
    Xe: (2, 3) array of triangle vertex coordinates, one column per vertex
  
  Returns:
    J: (2, 2) affine transformation matrix for reference triangle
  """
  return np.column_stack((Xe[:, 1] - Xe[:, 0], Xe[:, 2] - Xe[:, 0]))

def map_quad_grid(Rq, Sq, tris):
  """
  Maps quadrature points on the reference triangle
  to the triangles contained in tris
  
  Parameters:
    Rq, Sq: quadrature points R,S on T_ref
    tris: delaunay triangulation

  Returns:
    XX,YY: stacked grid of mapped quad nodes on triangulation
  """
  Nrs = len(Rq)
  nTri = len(tris.simplices)
  XX = np.zeros((Nrs * nTri,1))
  YY = np.zeros((Nrs * nTri,1))
  RS = np.stack([Rq,Sq], axis=1)
  for iTri in range(nTri):
    tri_ind = tris.simplices[iTri]
    Pts_Ti = tris.points[tri_ind].T
    J = incidence_matrix(Pts_Ti)
    XYe = (J @ RS.T + Pts_Ti[:,0].reshape(2,1)).T
    XX[iTri*Nrs:Nrs*(iTri+1)] = XYe[:,0].reshape(Nrs,1)
    YY[iTri*Nrs:Nrs*(iTri+1)] = XYe[:,1].reshape(Nrs,1)
  return XX, YY


# for using circles to model obstacle hits (not great)
def estimate_obstacle_radius(cone_angle_deg, max_range_cm, n_sensors):
  delta_theta = np.radians(cone_angle_deg) / (n_sensors - 1)
  return max_range_cm * np.sin(delta_theta / 2)


# sigma based on empirical data for sensor uncertainty at distance
def sigma_rule_of_thumb(distance_cm, scale = 1.0):
  """
  Returns the standard deviation (in cm) using the 1% rule of thumb.
  """
  return scale * np.maximum(0.3, 0.01 * distance_cm)

# gaussian splatter to model 
def rho_splatter(grid_x, grid_y, hit_points, sensor_points, eps=0.99):
  """
  Evaluates the obstacle field ρ(x, y) from Gaussian splatter at hit points.

  Parameters:
    grid_x, grid_y : 2D meshgrid arrays (same shape) in world coordinates (cm)
    hit_points     : list of (x, y) coordinates of sensor hits
    sensor_points  : corresponding list of (x_s, y_s) sensor origins
    eps            : peak amplitude of each Gaussian (default: 0.99)

  Returns:
    rho : 2D array with obstacle field values ρ(x, y) ∈ [0, 1]
  """
  rho = np.zeros_like(grid_x)
  scale = 0.1
  min_sigma = 0.02
  for (xi, yi), (xs, ys) in zip(hit_points, sensor_points):
    d = np.hypot(xi - xs, yi - ys)  # Euclidean distance (in cm)
    sigma = min_sigma + scale * d
    exponent = -((grid_x - xi)**2 + (grid_y - yi)**2) / (2 * sigma**2)
    rho += eps * np.exp(exponent)
  return np.clip(rho, 0, 1)
