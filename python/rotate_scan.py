#!/usr/bin/env python3
# rotate_scan.py — rotate by a target angle (radians) and plot sonar hit points
# Requires: AlphaBot2.py in PYTHONPATH, matplotlib installed
# Indentation = two spaces; angles in radians

import time
import math
import csv
import os
import RPi.GPIO as GPIO
import numpy as np
from AlphaBot2 import AlphaBot2
from scipy.spatial import Delaunay

# ---- USER SETTINGS ----------------------------------------------------------

# Ultrasonic pins (AlphaBot2 slot you just confirmed)
TRIG = 22
ECHO = 27

# Rotation calibration:
# Measure how long (seconds) it takes your bot to turn 90° (pi/2 rad) at a given speed.
# sec_per_rad = seconds_for_90deg / (math.pi / 2)
SECONDS_FOR_360_AT_SPEED = 1.1 # <-- put your measured value here
SEC_PER_RAD = SECONDS_FOR_360_AT_SPEED / (2.0 * math.pi)

# Motor speed (PWM duty, ~ -100..100). Keep modest to reduce slip.
ROTATE_SPEED = 20

# Scan settings
SAMPLE_PERIOD = 0.01  # seconds between sonar shots
ECHO_TIMEOUT = 0.03   # seconds (30 ms -> ~5 m max)
MIN_CM, MAX_CM = 2.0, 400.0

# Output
OUT_CSV = "scan_points.csv"
OUT_PNG = "scan_points.png"

# Direction convention:
# target_angle_rad > 0  -> clockwise turn (left wheel forward, right wheel backward)
# target_angle_rad < 0  -> counter-clockwise
# You can change this mapping if your bot spins the opposite way.

# ---- ULTRASONIC -------------------------------------------------------------

def setup_ultrasonic():
  GPIO.setmode(GPIO.BCM)
  GPIO.setup(TRIG, GPIO.OUT, initial=GPIO.LOW)
  GPIO.setup(ECHO, GPIO.IN, pull_up_down=GPIO.PUD_DOWN)
  time.sleep(0.05)  # settle

def read_distance_cm(timeout=ECHO_TIMEOUT):
  # 10 µs trigger pulse
  GPIO.output(TRIG, GPIO.HIGH)
  time.sleep(10e-6)
  GPIO.output(TRIG, GPIO.LOW)

  t0 = time.time()
  # wait for rising edge
  while GPIO.input(ECHO) == 0:
    if time.time() - t0 > timeout:
      return None
  ts = time.time()

  # wait for falling edge
  while GPIO.input(ECHO) == 1:
    if time.time() - ts > timeout:
      return None
  te = time.time()

  dur = te - ts
  dist_cm = (dur * 34300.0) / 2.0
  if dist_cm < MIN_CM or dist_cm > MAX_CM:
    return None
  return dist_cm

# ---- ROTATION & SCAN --------------------------------------------------------

def rotate_and_scan(target_angle_rad):
  bot = AlphaBot2()
  setup_ultrasonic()

  # choose direction
  if target_angle_rad > 0:
    # clockwise: left forward, right backward
    bot.setMotor(ROTATE_SPEED, -ROTATE_SPEED)
    direction = 1.0
  else:
    bot.setMotor(-ROTATE_SPEED, ROTATE_SPEED)
    direction = -1.0

  start = time.time()
  last_sample = 0.0
  points = []   # (theta_rad, dist_cm)

  try:
    while True:
      t = time.time() - start
      theta = direction * (t / SEC_PER_RAD)  # estimated rotation since start (rad)

      # sample sonar on schedule
      if t - last_sample >= SAMPLE_PERIOD:
        d = read_distance_cm()
        points.append((theta, d))
        last_sample = t

      # stop condition
      if abs(theta) >= abs(target_angle_rad):
        break
      # tiny sleep to reduce CPU
      time.sleep(0.002)
  finally:
    bot.stop()
    # --- return to original orientation ---
    elapsed = time.time() - start
    return_dir = -direction
    bot.setMotor(return_dir * ROTATE_SPEED, -return_dir * ROTATE_SPEED)
    time.sleep(elapsed)
    bot.stop()
    GPIO.cleanup()

  return points

# ---- SAVE & PLOT ------------------------------------------------------------

def save_points(points, csv_path=OUT_CSV):
  with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["theta_rad", "distance_cm", "x_m", "y_m"])
    for theta, d in points:
      if d is None:
        w.writerow([theta, "", "", ""])
      else:
        r_m = d / 100.0
        x = r_m * math.cos(theta)
        y = r_m * math.sin(theta)
        w.writerow([theta, d, x, y])
def manhattan_rectangle_from_points(xs_cm, ys_cm, buffer_cm=5.0):
  """
  Smallest axis-aligned rectangle (in the robot's local frame) that contains all points,
  expanded by 'buffer_cm' on all sides. Returns corners CCW (4,2) in cm.
  """
  if xs_cm.size == 0 or ys_cm.size == 0:
    x_min, x_max = -buffer_cm, buffer_cm
    y_min, y_max = -buffer_cm, buffer_cm
  else:
    x_min = float(np.min(xs_cm) - buffer_cm)
    x_max = float(np.max(xs_cm) + buffer_cm)
    y_min = float(np.min(ys_cm) - buffer_cm)
    y_max = float(np.max(ys_cm) + buffer_cm)

  return np.array([
    [x_min, y_min],  # bottom-left
    [x_min, y_max],  # top-left
    [x_max, y_max],  # top-right
    [x_max, y_min],  # bottom-right
  ], dtype=float)

def rect_bounds_from_corners(rect_local_cm):
  """
  Given CCW rectangle corners in local frame, return (x_min, x_max, y_min, y_max).
  """
  xs = rect_local_cm[:, 0]
  ys = rect_local_cm[:, 1]
  return float(xs.min()), float(xs.max()), float(ys.min()), float(ys.max())

def local_to_world(points_local_cm, robot_pose):
  """
  Map Nx2 points from robot local frame (cm) to world (cm).
  robot_pose = (x0_cm, y0_cm, theta_rad).
  """
  x0, y0, th = robot_pose
  R = np.array([[math.cos(th), -math.sin(th)],
                [math.sin(th),  math.cos(th)]], dtype=float)
  return (R @ points_local_cm.T).T + np.array([x0, y0], dtype=float)

def world_to_local(points_world_cm, robot_pose):
  """
  Map Nx2 points from world (cm) to robot local frame (cm).
  """
  x0, y0, th = robot_pose
  Rt = np.array([[ math.cos(th),  math.sin(th)],
                 [-math.sin(th),  math.cos(th)]], dtype=float)  # R^T
  return (points_world_cm - np.array([x0, y0], dtype=float)) @ Rt.T

def affine_world_to_canonical_rect(rect_world_cm, robot_pose, eps=1e-9):
  """
  Build an affine map [u,v]^T = A @ [Xw,Yw]^T + b that sends the *world-space*
  rectangle (image of a local axis-aligned rectangle under the pose) to
  the canonical domain [0,1] x [-1,1].

  Steps:
    world -> local (rigid)
    local (x in [x_min,x_max], y in [y_min,y_max]) -> canonical:
      u = (x - x_min) / (x_max - x_min)
      v = 2 * (y - y_min) / (y_max - y_min) - 1

  Returns:
    A (2x2), b (2,), and the local bounds (x_min, x_max, y_min, y_max)
  """
  # Recover local rectangle bounds by inverting pose on the given world corners
  rect_local = world_to_local(rect_world_cm, robot_pose)
  x_min, x_max, y_min, y_max = rect_bounds_from_corners(rect_local)

  Wx = max(x_max - x_min, eps)  # width  in local-x
  Hy = max(y_max - y_min, eps)  # height in local-y

  # local -> canonical
  # u = (x - x_min)/Wx
  # v = 2*(y - y_min)/Hy - 1
  A_local_to_can = np.array([[1.0/Wx,         0.0],
                             [0.0,     2.0/Hy]], dtype=float)
  b_local_to_can = np.array([-x_min/Wx, -2.0*y_min/Hy - 1.0], dtype=float)

  # world -> local
  x0, y0, th = robot_pose
  Rt = np.array([[ math.cos(th),  math.sin(th)],
                 [-math.sin(th),  math.cos(th)]], dtype=float)
  A_world_to_local = Rt
  b_world_to_local = -Rt @ np.array([x0, y0], dtype=float)

  # Compose: [u,v] = A_l2c @ (A_w2l @ [Xw,Yw] + b_w2l) + b_l2c
  A = A_local_to_can @ A_world_to_local
  b = A_local_to_can @ b_world_to_local + b_local_to_can
  return A, b, (x_min, x_max, y_min, y_max)

def canonical_to_world_rect(grid_can, robot_pose, local_bounds):
  """
  Inverse map: canonical -> world for an arbitrary axis-aligned local rectangle.
  Inputs:
    grid_can: (M,2) with columns (u in [0,1], v in [-1,1])
    robot_pose: (x0_cm, y0_cm, theta_rad)
    local_bounds: (x_min, x_max, y_min, y_max) in *local* cm
  Output:
    grid_world: (M,2) world cm
  """
  x_min, x_max, y_min, y_max = local_bounds
  Wx = x_max - x_min
  Hy = y_max - y_min

  u = grid_can[:, 0]
  v = grid_can[:, 1]

  # canonical -> local
  x_loc = x_min + u * Wx
  y_loc = y_min + 0.5 * (v + 1.0) * Hy
  pts_local = np.column_stack([x_loc, y_loc])

  # local -> world
  return local_to_world(pts_local, robot_pose)

def canonical_grid_and_triangulation(nu, nv=None):
  """
  Canonical domain grid: u in [0,1] with (nu+1) nodes,
  v in [-1,1] with (nv+1) nodes (defaults to 2*nu for aspect).
  Returns (grid_can (M,2), tri) with Delaunay over canonical nodes.
  """
  if nv is None:
    nv = nu
  u_lin = np.linspace(0.0, 1.0, nu + 1)
  v_lin = np.linspace(-1.0, 1.0, nv + 1)
  UU, VV = np.meshgrid(u_lin, v_lin)
  grid_can = np.column_stack([UU.ravel(), VV.ravel()])
  tri = Delaunay(grid_can)
  return grid_can, tri


def polar_hits_to_xy_cm(points):
  """Convert (theta_rad, dist_cm) -> arrays x_cm, y_cm (skip Nones)."""
  xs, ys = [], []
  for theta, d in points:
    if d is None: 
      continue
    r = d  # already in cm
    xs.append(r * math.cos(theta))
    ys.append(r * math.sin(theta))
  return np.asarray(xs), np.asarray(ys)

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


def plot_points(points, png_path=OUT_PNG):
  import matplotlib
  matplotlib.use("Agg")  # headless
  import matplotlib.pyplot as plt

  Rq = np.loadtxt("xtri_N496_n30_M1378_m51.txt");
  Sq = np.loadtxt("ytri_N496_n30_M1378_m51.txt");
  Wq = np.loadtxt("wtri_N496_n30_M1378_m51.txt"); 

  xs, ys = polar_hits_to_xy_cm(points)

  # robot pose is 0 for now
  robot_pose = (0., 0., 0.) 
  robot_x, robot_y, robot_theta = robot_pose
  # generate robot sweep box
  rect_local = manhattan_rectangle_from_points(xs, ys)
  R = np.array([
    [np.cos(robot_theta), -np.sin(robot_theta)],
    [np.sin(robot_theta),  np.cos(robot_theta)]
  ])
  rect_world = local_to_world(rect_local, robot_pose)
  fig, ax = plt.subplots(1,2)
  ax[0].plot(rect_world[:, 0], rect_world[:, 1], 'b-', label="Robot")
  ax[0].plot((rect_world[-1, 0], rect_world[0,0]),\
           (rect_world[-1, 1], rect_world[0,1]), 'b-', label="Robot")
  ax[0].scatter(xs, ys, s=10)
  ax[0].scatter([0], [0], marker="x")  # robot origin
  ax[0].set_title("Ultrasonic hits during rotation")
  ax[0].set_xlabel("x (m)")
  ax[0].set_ylabel("y (m)")
  ax[0].grid(True)


  A_w2c, b_w2c, local_bounds = affine_world_to_canonical_rect(rect_world, robot_pose)
  # Now any world point Pw maps as: [u,v]^T = A_w2c @ Pw + b_w2c 
  
  # canonical grid and triangulation
  grid_can, tris =  canonical_grid_and_triangulation(4)
  # quad points mapped to each triangle in sector
  XX, YY = map_quad_grid(Rq, Sq, tris)

  ax[1].scatter(XX,YY,s=10)
  ax[1].grid(True)
  ax[1].set_title("Quad points in canonical domain")
  ax[1].set_xlabel("x")
  ax[1].set_xlabel("y")
  fig.tight_layout()
  fig.savefig(png_path, dpi=160)


# ---- MAIN -------------------------------------------------------------------

if __name__ == "__main__":
  import argparse
  ap = argparse.ArgumentParser(description="Rotate AlphaBot2 by angle (rad) and plot HC-SR04 hits.")
  ap.add_argument("angle_rad", type=float, help="Target rotation (radians). Positive = clockwise.")
  args = ap.parse_args()

  pts = rotate_and_scan(args.angle_rad)
  save_points(pts, OUT_CSV)
  plot_points(pts, OUT_PNG)

  print(f"Saved {len(pts)} samples to {OUT_CSV} and plot to {OUT_PNG}.")
