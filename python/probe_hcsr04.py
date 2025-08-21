# probe_hcsr04.py
import time, sys
import RPi.GPIO as GPIO

# Safe set of GPIOs to try (BCM numbering), skipping 2/3 to avoid I2C conflicts.
outs = [4, 17, 27, 22, 23, 24, 5, 6, 12, 13, 16, 18, 19, 20, 21, 26]
ins  = [4, 17, 27, 22, 23, 24, 5, 6, 12, 13, 16, 18, 19, 20, 21, 26]

def try_pair(TRIG, ECHO, tries=3, timeout=0.03):
  GPIO.setmode(GPIO.BCM)
  # Start with everything low and echo pulled down
  GPIO.setup(TRIG, GPIO.OUT, initial=GPIO.LOW)
  GPIO.setup(ECHO, GPIO.IN, pull_up_down=GPIO.PUD_DOWN)
  time.sleep(0.05)
  vals = []
  for _ in range(tries):
    # 10 us TRIG pulse
    GPIO.output(TRIG, GPIO.HIGH)
    time.sleep(10e-6)
    GPIO.output(TRIG, GPIO.LOW)

    # wait for rising edge
    t0 = time.time()
    while GPIO.input(ECHO) == 0:
      if time.time() - t0 > timeout:
        vals.append(None)
        break
    else:
      ts = time.time()
      # wait for falling edge
      while GPIO.input(ECHO) == 1:
        if time.time() - ts > timeout:
          vals.append(None)
          break
      else:
        te = time.time()
        dist = (te - ts) * 34300 / 2
        vals.append(dist)
    time.sleep(0.06)
  GPIO.cleanup()
  # Count successes
  good = [v for v in vals if v is not None and 1.0 <= v <= 400.0]
  return vals, len(good)

tested = 0
for trig in outs:
  for echo in ins:
    if trig == echo:  # avoid same pin
      continue
    vals, ok = try_pair(trig, echo)
    tested += 1
    if ok >= 2:
      print(f"LIKELY MATCH: TRIG={trig} ECHO={echo} -> {['timeout' if v is None else round(v,1) for v in vals]}")
      sys.exit(0)
print(f"Scanned {tested} pairs; no match found.")
