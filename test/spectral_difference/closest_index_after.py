import numpy as np
import spectral_difference as sd

filtered=np.array([2,4,6,8])
find_closest=np.array([0,3,7,9])
left_closest=sd.closest_index_after(filtered,find_closest,reverse=True)

pairs=sd.up_down_match(filtered,left_closest)
print('attacks')
print(filtered)
print('left closest')
print(left_closest)
print('pairs')
print(pairs)
