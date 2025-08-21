import time
import RPi.GPIO as GPIO

candidates = [
  (27, 22),
  (23, 24),
  (17, 18),
  (5, 6),
  (12, 13),  # uncommon, but some boards repurpose
  (16, 20),  # ditto
]

def measure(trig, echo, tries=4):
  GPIO.setmode(GPIO.BCM)
  GPIO.setup(trig, GPIO.OUT, initial=GPIO.LOW)
  GPIO.setup(echo, GPIO.IN, pull_up_down=GPIO.PUD_DOWN)
  time.sleep(0.05)
  ok = 0
  vals = []
  for _ in range(tries):
    # 10us pulse
    GPIO.output(trig, GPIO.HIGH)
    time.sleep(10e-6)
    GPIO.output(trig, GPIO.LOW)

    t0 = time.time()
    # wait for rising
    while GPIO.input(echo) == 0:
      if time.time() - t0 > 0.03:
        vals.append(None)
        break
    else:
      ts = time.time()
      # wait for falling
      while GPIO.input(echo) == 1:
        if time.time() - ts > 0.03:
          vals.append(None)
          break
      else:
        te = time.time()
        dist = (te - ts) * 34300 / 2
        vals.append(dist)
        ok += 1
    time.sleep(0.08)
  GPIO.cleanup()
  return ok, vals

for trig, echo in candidates:
  ok, vals = measure(trig, echo)
  print(f"TRIG {trig} / ECHO {echo} -> {vals}")
