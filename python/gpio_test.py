import RPi.GPIO as GPIO
import time

# Pin setup (BCM mode)
TRIG_1 = 23  # GPIO23, physical pin 16
ECHO_1 = 24  # GPIO24, physical pin 18

TRIG_2 = 5   # GPIO5, physical pin 29
ECHO_2 = 6   # GPIO6, physical pin 31

GPIO.setmode(GPIO.BCM)
GPIO.setup(TRIG_1, GPIO.OUT)
GPIO.setup(ECHO_1, GPIO.IN)
GPIO.setup(TRIG_2, GPIO.OUT)
GPIO.setup(ECHO_2, GPIO.IN)

def measure_distance(trig_pin, echo_pin, sensor_name="", max_distance=400.0):
  # Send 10µs trigger pulse
  GPIO.output(trig_pin, False)
  time.sleep(0.05)
  GPIO.output(trig_pin, True)
  time.sleep(0.00001)
  GPIO.output(trig_pin, False)

  # Wait for echo to start
  timeout = time.time() + 0.02
  while GPIO.input(echo_pin) == 0:
    pulse_start = time.time()
    if pulse_start > timeout:
      print(f"{sensor_name} Timeout: No echo start → returning max range")
      return max_distance

  # Wait for echo to end
  timeout = time.time() + 0.02
  while GPIO.input(echo_pin) == 1:
    pulse_end = time.time()
    if pulse_end > timeout:
      print(f"{sensor_name} Timeout: No echo end → returning max range")
      return max_distance

  duration = pulse_end - pulse_start
  distance = round(duration * 17150, 2)
  print(f"{sensor_name} Distance: {distance} cm")
  return distance

try:
  while True:
    measure_distance(TRIG_1, ECHO_1, "Sensor 1")
    time.sleep(0.1)  # Give time between sensors
    measure_distance(TRIG_2, ECHO_2, "Sensor 2")
    time.sleep(1)

except KeyboardInterrupt:
  print("\nExiting...")
  GPIO.cleanup()
