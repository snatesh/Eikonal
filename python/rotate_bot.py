import time
from AlphaBot2 import AlphaBot2

bot = AlphaBot2()

try:
    # Rotate in place: left motor forward, right motor backward
    bot.setMotor(20, -20)   # values = PWM duty cycle, adjust speed as needed
    time.sleep(1.25)         # adjust this duration for angle turned
    bot.stop()
finally:
    bot.stop()
