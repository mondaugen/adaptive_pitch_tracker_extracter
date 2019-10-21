import numpy as np
from scipy import signal
import subprocess

N=1000000
x=np.random.standard_normal(N).astype('float32')
x.tofile('/tmp/in.f32')
# do computation with C implementation
subprocess.run('src/test/bin/iir_lattice_filter')
# read in result
y_c=np.fromfile('/tmp/out.f32',dtype='float32')
a=[1,-0.8,0.64,-0.512]
b=[.328]
y=signal.lfilter(b,a,x)
# trim result
y_c=y_c[:len(y)]
print(np.any(np.isnan(y_c)))
print(np.sum(np.abs(y-y_c))/N)
print(y[:10])
print(y_c[:10])
