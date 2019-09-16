# create /tmp/guit.f64 via
# sox sounds/guitar.wav -c 1 -r 16k /tmp/guit.f64

import numpy as np
import classic_puckette_timestretch
import common

x=np.fromfile('/tmp/guit.f64',dtype='float64')

H=128
W=1024
N=len(x)
# simple 2x slower
t_2x=np.concatenate(([0],np.arange(H,N-W,H/2))).astype('int')

y_2x=classic_puckette_timestretch.time_stretch_arb_times(x,t_2x,H,W)

y_2x=common.normalize(y_2x)
y_2x.tofile('/tmp/guit_2x.f64')

# random
t_rand=[0]
t_ = H
while t_ <= N-W:
    t_rand.append(t_)
    t_+=np.random.randint(2*H)
t_rand=np.array(t_rand)
t_rand=np.random.randint(H,N-W,10000)
t_rand=np.concatenate(([0],t_rand))

y_rand=classic_puckette_timestretch.time_stretch_arb_times(x,t_rand,H,W)
y_rand.tofile('/tmp/guit_rand.f64')
