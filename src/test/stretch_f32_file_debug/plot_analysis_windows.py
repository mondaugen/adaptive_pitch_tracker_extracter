import numpy as np
import matplotlib.pyplot as plt

an=np.fromfile("/tmp/an.f32",dtype='float32')
sy=np.fromfile("/tmp/sy.f32",dtype='float32')

fig,ax=plt.subplots(2,1)
ax[0].set_title('analysis')
ax[0].plot(np.arange(len(an)),an)
ax[1].set_title('synthesis')
ax[1].plot(np.arange(len(sy)),sy)

plt.show()

