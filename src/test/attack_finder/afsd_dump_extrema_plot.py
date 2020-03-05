import numpy as np
import matplotlib.pyplot as plt

sd_mins=np.fromfile('/tmp/sd_mins.u32',dtype='uint32')
sd_maxs=np.fromfile('/tmp/sd_maxs.u32',dtype='uint32')
spec_diff=np.fromfile('/tmp/spec_diff.f32',dtype='float32')
n=np.arange(len(spec_diff))
plt.plot(n,-1*spec_diff)
plt.plot(sd_mins,spec_diff[sd_mins],'.',label='mins')
plt.plot(sd_maxs,spec_diff[sd_maxs],'.',label='maxs')
plt.legend()

plt.show()
