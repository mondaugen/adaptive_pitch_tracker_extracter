import bl_square
import matplotlib.pyplot as plt
import numpy as np

bls=bl_square.bl_square_synth(D=0.9,normalize=True)

fig,axs=plt.subplots(2,1)

axs[0].plot(np.arange(bls.N),bls.tabs[0,:])
axs[1].plot(np.arange(bls.N),bls.tabs[len(bls.tabs)-3,:])
print(bls.n_tabs)

plt.show()
