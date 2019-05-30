import ptracker

import matplotlib.pyplot as plt

fig,ax=plt.subplots(1,1)

cf = ptracker.pitch_check_comb(50/44100,0.999)
cf.plot(ax)

plt.show()
