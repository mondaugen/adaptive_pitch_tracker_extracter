import numpy as np

N=2048
P=250
scale=2**16
i=scale//16
x=np.arange(N).astype('float32')
y=np.cos(2*np.pi*x/P).astype('float32')
xi=np.arange(scale,scale*(N-2),i,dtype='uint32')
xr=xi>>16

x.tofile('/tmp/x.f32')
y.tofile('/tmp/y.f32')
xi.tofile('/tmp/xi.u16q16')
xr.tofile('/tmp/xr.u32')
