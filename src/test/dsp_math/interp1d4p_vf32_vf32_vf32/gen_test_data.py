import numpy as np

N=2048
P=250
i=1/16
x=np.arange(N).astype('float32')
y=np.cos(2*np.pi*x/P).astype('float32')
xi=np.arange(1,N-2,i).astype('float32')
xr=np.round(xi).astype('uint32')

x.tofile('/tmp/x.f32')
y.tofile('/tmp/y.f32')
xi.tofile('/tmp/xi.f32')
xr.tofile('/tmp/xr.u32')
