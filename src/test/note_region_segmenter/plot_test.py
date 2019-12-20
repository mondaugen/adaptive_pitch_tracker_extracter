import numpy as np
import matplotlib.pyplot as plt
import common

B=256

active=np.fromfile('/tmp/active.f32',dtype='float32')
start=np.fromfile('/tmp/start.f32',dtype='float32')
end=np.fromfile('/tmp/end.f32',dtype='float32')
N=len(active)
n=np.arange(N)

common.logic_plot(n,start)
common.logic_plot(n,end+1)
common.logic_plot(n,active+2)

regions=np.fromfile('/tmp/regions.txt',sep=' ').reshape(-1,2)
for reg in regions:
    common.region_plot(reg[0],reg[1],height=0.5,level=3,color='k')
for n in range(0,N,B):
    common.region_plot(n,n+B,height=0.5,level=3.5,color='k')

plt.show()

