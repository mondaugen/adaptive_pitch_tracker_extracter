import numpy as np
import matplotlib.pyplot as plt
import spectral_difference


filtered=np.fromfile('/tmp/closest_index_after_test_filtered_forward.u32',dtype='uint32')
find_closest=np.fromfile('/tmp/closest_index_after_test_find_closest_forward.u32',dtype='uint32')
forward=np.fromfile('/tmp/closest_index_after_test_closest_forward.u32',dtype='uint32')
reverse=np.fromfile('/tmp/closest_index_after_test_closest_reverse.u32',dtype='uint32')

forward_truth=spectral_difference.closest_index_after(filtered,find_closest,reverse=False)
reverse_truth=spectral_difference.closest_index_after(filtered,find_closest,reverse=True)

fig,axs=plt.subplots(3,2)

axs[0,0].plot(filtered,np.ones_like(filtered),'.',label='filtered')
axs[0,0].plot(find_closest,np.ones_like(find_closest),'.',label='find_closest')
axs[0,0].legend()
axs[1,0].plot(forward,np.ones_like(forward),'.',label='forward')
axs[1,0].set_title('forward')
axs[2,0].plot(reverse,np.ones_like(reverse),'.',label='reverse')
axs[2,0].set_title('reverse')

axs[1,1].plot(forward_truth,np.ones_like(forward_truth),'.',label='forward_truth')
axs[1,1].set_title('forward_truth')
axs[2,1].plot(reverse_truth,np.ones_like(reverse_truth),'.',label='reverse_truth')
axs[2,1].set_title('reverse_truth')

plt.show()
