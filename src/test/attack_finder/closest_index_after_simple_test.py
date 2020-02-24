import numpy as np

filtered=np.array([10,20,30,40],dtype='uint32')
find_closest=np.array([15,31,34],dtype='uint32')
filtered.tofile('/tmp/closest_index_after_test_filtered_forward.u32')
filtered.tofile('/tmp/closest_index_after_test_filtered_reverse.u32')
find_closest.tofile('/tmp/closest_index_after_test_find_closest_forward.u32')
find_closest.tofile('/tmp/closest_index_after_test_find_closest_reverse.u32')
