import numpy as np
from scipy import signal
import subprocess
import filters

N=1000000
x=np.random.standard_normal(N).astype('float32')
x.tofile('/tmp/in.f32')
b=np.array([0.2759,0.5121,0.5121,0.2759])
a=np.array([1,-0.0010,.6546,-.0775])
r=filters.a_to_r(a)
print('r')
print(r)
A=filters.r_to_A(r)
print('A')
print(A)
c=filters.b_A_to_c(b,A)
print('c')
print(c)
r=r.astype('float32')
c=c.astype('float32')
r.tofile('/tmp/r.f32')
c=c.tofile('/tmp/c.f32')
# do computation with C implementation
env=dict(P='%d'%(len(r),))
subprocess.run('src/test/bin/lattice_filter_proc',env=env)
# read in result
y_c=np.fromfile('/tmp/out.f32',dtype='float32')
y=signal.lfilter(b,a,x)
# trim result
y_c=y_c[:len(y)]
print("Any result of C implementation NaN?",np.any(np.isnan(y_c)))
print("Average absolute difference between C and Python implementation:",
np.sum(np.abs(y-y_c))/N)
print("Python implementation")
print(y[:10])
print("C implementation")
print(y_c[:10])
