import numpy as np
import matplotlib.pyplot as plt
from os import environ

W=int(environ['W'])

W=W//2+1
# because we use real data the length is W//2+1

z_input0=np.fromfile("/tmp/z_input0.z64",dtype="complex64")
z_inputH=np.fromfile("/tmp/z_inputH.z64",dtype="complex64")

zi0=np.reshape(z_input0,(-1,W))
ziH=np.reshape(z_inputH,(-1,W))

#zi0b=np.isnan(zi0).astype('int')
#ziHb=np.isnan(ziH).astype('int')
zi0b=(complex(0) == zi0).astype('int')
ziHb=(complex(0) == ziH).astype('int')

zi0m=20*np.log10(np.abs(zi0))
ziHm=20*np.log10(np.abs(ziH))

fig,ax=plt.subplots(4,1)

ax[0].imshow(zi0b,aspect='auto',origin='lower')
ax[0].set_title('z_input0')
ax[1].imshow(ziHb,aspect='auto',origin='lower')
ax[1].set_title('z_inputH')
ax[2].imshow(zi0m,aspect='auto',origin='lower')
ax[3].imshow(ziHm,aspect='auto',origin='lower')

plt.show()
