import numpy as np
import matplotlib.pyplot as plt
import common

B=256

adsr_active=np.fromfile('/tmp/adsr_active.f32',dtype='float32')
start=np.fromfile('/tmp/start.f32',dtype='float32')
end=np.fromfile('/tmp/end.f32',dtype='float32')
adsr_envelope=np.fromfile('/tmp/adsr_envelope.f32',dtype='float32')
P=len(adsr_envelope)
N=len(adsr_active)
n=np.arange(N)

common.logic_plot(n,start,label='start')
common.logic_plot(n,end+1,label='end')
common.logic_plot(n,adsr_active+2,label='adsr_active')

regions=np.fromfile('/tmp/regions.txt',sep=' ').reshape(-1,2)
for reg in regions:
    common.region_plot(reg[0],reg[1],height=0.5,level=3,color='k')
for m in range(0,P,B):
    common.region_plot(m,m+B,height=0.5,level=3.5,color='k')
p=np.arange(P)
plt.plot(p,adsr_envelope+4)
plt.legend()

plt.show()

