#!/usr/bin/env python3
# rotate_scan.py — rotate by a target angle (radians) and plot sonar hit points
# Requires: AlphaBot2.py in PYTHONPATH, matplotlib installed
# Indentation = two spaces; angles in radians

import time
import math
import csv
import os
import RPi.GPIO as GPIO
from AlphaBot2 import AlphaBot2

# ---- USER SETTINGS ----------------------------------------------------------

# Ultrasonic pins (AlphaBot2 slot you just confirmed)
TRIG = 22
ECHO = 27

# Rotation calibration:
# Measure how long (seconds) it takes your bot to turn 90° (pi/2 rad) at a given speed.
# sec_per_rad = seconds_for_90deg / (math.pi / 2)
SECONDS_FOR_90_AT_SPEED = 0.375 # <-- put your measured value here
SEC_PER_RAD = SECONDS_FOR_90_AT_SPEED / (math.pi / 2)

# Motor speed (PWM duty, ~ -100..100). Keep modest to reduce slip.
ROTATE_SPEED = 20

# Scan settings
SAMPLE_PERIOD = 0.05  # seconds between sonar shots
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

def plot_points(points, png_path=OUT_PNG):
  import matplotlib
  matplotlib.use("Agg")  # headless
  import matplotlib.pyplot as plt

  xs, ys = [], []
  for theta, d in points:
    if d is None:
      continue
    r_m = d / 100.0
    xs.append(r_m * math.cos(theta))
    ys.append(r_m * math.sin(theta))

  plt.figure(figsize=(5, 5))
  plt.scatter(xs, ys, s=10)
  plt.scatter([0], [0], marker="x")  # robot origin
  plt.gca().set_aspect("equal", adjustable="box")
  plt.title("Ultrasonic hits during rotation")
  plt.xlabel("x (m)")
  plt.ylabel("y (m)")
  plt.grid(True)
  plt.tight_layout()
  plt.savefig(png_path, dpi=160)

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
  print("Tip: tune SEC_PER_RAD by timing a known turn (e.g., pi/2).")
