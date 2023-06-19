import math
import numpy as np

pi = math.pi

degree = 30.0

dx = math.cos(degree*pi/180.0)
dy = math.sin(degree*pi/180.0)

print(dx, dy)


xyzi = [1,0,0]
xyzi = [0,1,0]
xyzi = [0,0,1]

ri = np.sqrt(np.sum(np.power(xyzi,2)))
theta  = math.acos(xyzi[2]/ri)
phi    = math.atan2(xyzi[1], xyzi[0])

print(ri, theta, phi)
