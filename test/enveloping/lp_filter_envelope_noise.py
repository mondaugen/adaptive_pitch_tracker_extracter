import numpy as np
from scipy import signal, interpolate
import filters
import common
import matplotlib.pyplot as plt

SR=16e3
N=100000
filter_order=4
x=np.random.standard_normal(N)
f_max=3000
f_sus=500

def filt_des(n,f):
    b,a=signal.butter(n,f,fs=1.)
    return (b,a)

filter_tab=filters.filter_interp_table(
        filt_des,
        filter_order,
        min_f=10/SR,
        out_of_bounds_values='lowpass',
        designs_per_octave=12,
        check_filter_stability=True)

filter_env=interpolate.interp1d(
    np.array([0,1000,10000,50000,N]),
    np.array([0,f_max,f_sus,f_sus,50]),
    fill_value=(0,0))(np.arange(N))/SR

R,C=filter_tab.lookup(filter_env)

y=filters.lattice_filter_proc(x,R,C)
print('beginning',y[:10])
print('end',y[-20:])
y=common.normalize(y)
y.tofile('/tmp/out.f32')

fig,axs=plt.subplots(3,1,sharex=True)
t=np.arange(N)/SR
axs[0].plot(t,y)
for n,r in enumerate(R):
    axs[1].plot(t,r,label='%d'%(n,))

axs[1].legend()

for n,c in enumerate(C):
    axs[2].plot(t,c,label='%d'%(n,))

axs[2].legend()

plt.show()

