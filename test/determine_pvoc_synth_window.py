import classic_puckette_timestretch as cpt
import numpy as np


H=8
W=32
x=np.random.standard_normal(1000)
a_times=np.arange(0,len(x)-W,H)
print(a_times)

y=cpt.time_stretch_arb_times(x,a_times,H,W,window_type='parzen')

print(y/x[:len(y)])
