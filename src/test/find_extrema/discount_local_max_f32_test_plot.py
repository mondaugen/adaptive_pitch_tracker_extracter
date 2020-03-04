import numpy as np
import matplotlib.pyplot as plt

x=np.fromfile('/tmp/noi.f32',dtype='float32')
thresh=np.fromfile('/tmp/noithresh.f32',dtype='float32')
maxs=np.fromfile('/tmp/noimaxs.u32',dtype='uint32')

print('n_thresh',len(thresh))
n_thresh=np.arange(len(thresh))
print('n_maxs',len(maxs))
n_maxs=np.arange(len(maxs))

plt.plot(n_thresh,thresh,label='thresh')
plt.plot(n_thresh,x,label='x')
plt.plot(maxs,x[maxs],'.',label='maxs')
plt.legend()

plt.show()
