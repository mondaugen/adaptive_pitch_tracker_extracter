import classic_puckette_timestretch as cpt
import numpy as np


H=5
W=32
x=np.random.standard_normal(1000)
a_times=np.arange(0,len(x)-W,H)
print('hop times')
print(a_times)

y=cpt.time_stretch_arb_times(x,a_times,H,W,
window_type='blackmanharris',
synth_window_type='hann')

print(y/x[:len(y)])
