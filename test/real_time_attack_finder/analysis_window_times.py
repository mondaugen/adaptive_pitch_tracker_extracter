# checking to see the analysis windows don't cross an attack

import numpy as np
import matplotlib.pyplot as plt
import attack_finder
import rtaf_common
import time_map_tstretch
import common

W=1024
H=256
# TODO: The ringbuffer overflows if M=W
# It seems with W//2 that it doesn't ignore the correct attack
M=W//2
rtaac=time_map_tstretch.real_time_attack_avoid_controller(M,W,H)
print(rtaac.R/H)
afsd=rtaf_common.afdf_rt_test(
        # hop size or block size
        H=H,
        # window size for spectral difference
        W=W,
        N=25600,
        attack_freq_limit=rtaac.R//H+2)

rh=-1
for h in range(0,afsd.N,afsd.H):
    afsd.attacks[h//afsd.H]=afsd.afsd(afsd.x[h:h+afsd.H])
    ai=-1
    if afsd.attacks[h//afsd.H] > 0:
        ai=0
    print(h)
    _,read_idx,reset=rtaac.write_read_cycle(afsd.x[h:h+afsd.H],ai)
    color = 'red' if reset else 'black'
    common.plot_arch(plt,read_idx,W+H,rh,color=color)
    rh-=1

h=np.arange(0,afsd.N,afsd.H)
n=np.arange(afsd.N)

plt.plot(n,np.abs(afsd.x),label='original')
plt.plot(h,afsd.attacks,'.',label='attacks')
plt.plot(h,afsd.thresh,label='thresh')
plt.plot(h,afsd.sd,label='spectral difference')
plt.legend()
plt.show()
