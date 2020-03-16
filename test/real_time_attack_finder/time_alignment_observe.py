# Checking the time alignment of the real time attack finder

import numpy as np
import matplotlib.pyplot as plt
import attack_finder

# hop size or block size
H=256
# window size for spectral difference
W=1024

N=100*H
n=np.arange(N)

x=np.zeros(N)
x[5000:15000:2000]=1
attacks=np.zeros(N//H)
thresh=np.zeros_like(attacks)
sd=np.zeros_like(attacks)

class record_thresh:
    def __init__(self):
        self.n=0
    def __call__(self,t):
        thresh[self.n]=t
        self.n += 1
rt=record_thresh()

class record_sd:
    def __init__(self):
        self.n=0
    def __call__(self,t):
        sd[self.n]=t
        self.n += 1
rsd=record_sd()

afsd=attack_finder.attacks_from_spectral_diff_rt(
W=W,
H=H,
record_thresh=rt,
record_sd=rsd,
lmax_filt_rate=500/H,
attack_freq_limit=int(np.ceil(3000/H)),
ng_th=-40)

for h in range(0,N,H):
    attacks[h//H]=afsd(x[h:h+H])

h=np.arange(0,N,H)

plt.plot(n,np.abs(x),label='original')
plt.plot(h,attacks,'.',label='attacks')
plt.plot(h,thresh,label='thresh')
plt.plot(h,sd,label='spectral difference')
plt.legend()
plt.show()
