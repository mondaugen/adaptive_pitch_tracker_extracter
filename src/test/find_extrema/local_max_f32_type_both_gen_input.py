import numpy as np

N=100000
n=np.arange(N)
x=np.cos(2*np.pi*n/1000)*np.cos(2*np.pi*n/50)
x.astype('float32').tofile("/tmp/local_max_f32_input.f32")
