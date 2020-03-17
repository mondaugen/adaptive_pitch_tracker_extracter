# Checking the time alignment of the real time attack finder

import numpy as np
import matplotlib.pyplot as plt
import attack_finder
import rtaf_common

H=256
W=1024
M=W
R=2*M+W+H
afsd=rtaf_common.afdf_rt_test(
H=H,
W=W,
a_dist=R+1,
last_a=150000,
attack_freq_limit=int(np.ceil((R+1)/H))
)
for h in range(0,afsd.N,afsd.H):
    afsd.attacks[h//afsd.H]=afsd.afsd(afsd.x[h:h+afsd.H])

h=np.arange(0,afsd.N,afsd.H)
n=np.arange(afsd.N)

plt.plot(n,np.abs(afsd.x),label='original')
plt.plot(h,afsd.attacks,'.',label='attacks')
plt.plot(h,afsd.thresh,label='thresh')
plt.plot(h,afsd.sd,label='spectral difference')
plt.legend()
plt.show()
