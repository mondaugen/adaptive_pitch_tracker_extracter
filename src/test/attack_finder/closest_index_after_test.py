import numpy as np

SEED=1234

N=10
MAX=100

filtered=np.random.random_integers(0,MAX,N).astype('uint32')
filtered.sort()
find_closest=np.random.random_integers(0,MAX,2*N).astype('uint32')
find_closest.sort()

filtered.tofile('/tmp/closest_index_after_test_filtered_forward.u32')
filtered.tofile('/tmp/closest_index_after_test_filtered_reverse.u32')
find_closest.tofile('/tmp/closest_index_after_test_find_closest_forward.u32')
find_closest.tofile('/tmp/closest_index_after_test_find_closest_reverse.u32')
