# test avoiding attacks in real time
import numpy as np
import time_map_tstretch
import matplotlib.pyplot as plt

# plot an arch showing where the analysis window is placed
def plot_arch(ax,s,l,h,b=0,color='black'):
    ax.plot([s,s],[b,b+h],color=color)
    ax.plot([s,s+l],[b+h,b+h],color=color)
    ax.plot([s+l,s+l],[b+h,b],color=color)

# margin on either side of attack
M = 1
# analysis window length
W = 16
# hop size
H = 4
# size of safe region
R = 2*M+W+H
# analysis hop size
Ha=1

# number of samples in test signal
N = 40*H
# dummy signal
x=np.arange(N)
x_attacks=np.zeros(N)
x_attacks[5] = 1
x_attacks[5+R+1] = 1
x_attacks[5+R+1+3*R] = 1

rtaac=time_map_tstretch.real_time_attack_avoid_controller(M,W,H)

read_times=np.zeros(N//Ha)
resets=np.zeros_like(read_times)

h=0
n_h=0
for n in range(0,N,H):
    attack_idx = time_map_tstretch.compute_attack_index(x_attacks[n:n+H])
    attack_idx = -1 if attack_idx is None else attack_idx
    rtaac.write(x[n:n+H],attack_idx)
    while n_h < (n+H):
        _,t,r=rtaac.read()
        read_times[h]=t
        resets[h]=float(r)
        h += 1
        n_h += Ha

plt.plot(np.arange(N),x_attacks,'.')
h=-1
for t,r in zip(read_times,resets):
    color = 'red' if r else 'black'
    plot_arch(plt,t,W+H,h,color=color)
    h-=1

plt.show()

