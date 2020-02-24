import numpy as np
import matplotlib.pyplot as plt
from matplotlib import markers
import spectral_difference


filtered=np.fromfile('/tmp/closest_index_after_test_filtered_forward.u32',dtype='uint32')
find_closest=np.fromfile('/tmp/closest_index_after_test_find_closest_forward.u32',dtype='uint32')
forward=np.fromfile('/tmp/closest_index_after_test_closest_forward.u32',dtype='uint32')
print('forward',forward)
reverse=np.fromfile('/tmp/closest_index_after_test_closest_reverse.u32',dtype='uint32')
print('reverse',reverse)

forward_truth=spectral_difference.closest_index_after(filtered,find_closest,reverse=False)
reverse_truth=spectral_difference.closest_index_after(filtered,find_closest,reverse=True)

fig,axs=plt.subplots(2,1)

for i in range(2):
    axs[i].plot(filtered,np.ones_like(filtered),marker=markers.TICKUP,label='filtered',ls='')
    axs[i].plot(find_closest,np.ones_like(find_closest),marker=markers.TICKDOWN,
    label='find_closest',ls='')

axs[0].plot(forward,np.ones_like(forward),marker=markers.CARETUP,ls='',label='estimated forward')
axs[0].plot(reverse,np.ones_like(reverse),marker=markers.CARETDOWN,ls='',label='estimated reverse')

axs[1].plot(forward_truth,np.ones_like(forward_truth),marker=markers.CARETUP,label='forward_truth',ls='')
axs[1].plot(reverse_truth,np.ones_like(reverse_truth),marker=markers.CARETDOWN,label='reverse_truth',ls='')

for i in range(2):
    axs[i].legend()

plt.show()
