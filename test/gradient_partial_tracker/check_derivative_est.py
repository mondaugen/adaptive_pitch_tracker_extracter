# Check the derivative estimates for each of the estimation styles
# (inner-product, inner-product with approximate ramp, recursive with
# approximate ramp).

import dft_hill_climbing as dhc
import recursive_dft
import numpy as np
from scipy import signal

N=10000
x=np.random.standard_normal(N)
v=0.01
Nw=4096
w=signal.get_window('hann',Nw)
# normalize
w/=np.sum(w)

# Direct method of computation
Xv_direct=dhc.dft_bin(x[:Nw]*w,v)
dvXv_direct= dhc.dft_bin_dv(x[:Nw]*w,v)
# Direct method but with exponential approximating the ramp
dvXv_approx=dhc.dft_bin_dv_approx(Nw,1e-4)
dvXv_exp=dvXv_approx(x[:Nw]*w,v)
# Recursive method
wp=np.array([0.5,-0.5]) # hann window
# normalize
wp/=(wp[0]*Nw)
dvXv_rec_computer=recursive_dft.rec_dv_dft_sumcos(Nw,wp,1,max_err=1e-6)
dvXv_rec=0
for n in range(Nw):
    Xv_rec,dvXv_rec=dvXv_rec_computer.update(x[n],np.array([v]))
    
print('dvXv_direct:',dvXv_direct)
print('dvXv_exp:',dvXv_exp)
print('dvXv_rec:',dvXv_rec)
print('dvXv_exp/dvXv_direct:',np.abs(dvXv_exp/dvXv_direct))
print('|dvXv_exp-dvXv_direct|:',np.abs(dvXv_exp-dvXv_direct))
print('dvXv_rec/dvXv_direct:',np.abs(dvXv_rec/dvXv_direct))
print('|dvXv_rec-dvXv_direct|:',np.abs(dvXv_rec-dvXv_direct))

print('Xv_direct:',Xv_direct)
print('Xv_rec:',Xv_rec)
