import numpy as np
from window_tools import signal_windower

sw=signal_windower(np.array([0,1,2,3,4,5,6,7,8,9]))

print(sw[0:5:2])
print(sw[-2:5:2])
print(sw[-4:5:2])
print(sw[-6:14:2])
print(sw[-6:14:1])
print(sw[0])
print(sw[1])
print(sw[-1])
print(sw[10])
