import time
from AlphaBot2 import AlphaBot2

def turn(bot, angle_deg, speed=20):
    # Calibrate your own factor (sec per degree at given speed)
    sec_per_deg = 1.5/360.0   # adjust from your measurements
    duration = abs(angle_deg) * sec_per_deg
    if angle_deg > 0:   # positive = clockwise
        bot.setMotor(speed, -speed)
    else:               # negative = counter-clockwise
        bot.setMotor(-speed, speed)
    time.sleep(duration)
    bot.stop()

# Example usage:
bot = AlphaBot2()
#turn(bot, 360)     # rotate ~90° clockwise
turn(bot, 360)    # rotate ~45° counter-clockwise
