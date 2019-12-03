import numpy as np
N=10
x=np.zeros((N,N,N),dtype='float32')
for i in range(N):
    for j in range(N):
        for k in range(N):
            x[i,j,k]=N*N*i+N*j+k
x.tofile('/tmp/multi_dims.f32')
