import numpy as np
from scipy import signal
import subprocess
import filters
from scipy import interpolate

f0=0.4
f1=0.0001
q0=300
q1=300
b0,a0=signal.iirpeak(f0,q0,fs=1)
r0,c0=filters.b_a_to_r_c(b0,a0)
b1,a1=signal.iirpeak(f1,q1,fs=1)
r1,c1=filters.b_a_to_r_c(b1,a1)
N=100000
def get_interp(a0,a1):
    return interpolate.interp1d(
        [0,N],
        np.concatenate((a0[:,None],a1[:,None]),axis=1))(np.arange(N))
r=get_interp(r0,r1)
c=get_interp(c0,c1)

x=np.random.standard_normal(N).astype('float32')
x.tofile('/tmp/in.f32')
r=r.astype('float32')
c=c.astype('float32')
r.T.tofile('/tmp/r.f32')
c.T.tofile('/tmp/c.f32')
env=dict(P='%d'%(len(r),))
subprocess.run('src/test/bin/lattice_filter_proc',env=env)
