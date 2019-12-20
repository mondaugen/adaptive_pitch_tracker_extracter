import numpy as np
import matplotlib.pyplot as plt
import common

active=np.fromfile('/tmp/active.f32',dtype='float32')
start=np.fromfile('/tmp/start.f32',dtype='float32')
end=np.fromfile('/tmp/end.f32',dtype='float32')
n=np.arange(len(active))

common.logic_plot(n,active)
common.logic_plot(n,active+1)
common.logic_plot(n,active+2)

regions=np.fromfile('/tmp/regions.txt',sep=' ').reshape(-1,2)
for reg in regions:
    common.region_plot(reg[0],reg[1],level=3,color='k')

plt.show()

