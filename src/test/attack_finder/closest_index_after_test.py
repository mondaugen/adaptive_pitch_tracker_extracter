import numpy as np
import random

SEED=1234

N=10
MAX=100

random.seed(SEED)
# values must be unique
filtered=np.array(random.sample(range(MAX),N),dtype='uint32')
filtered.sort()
find_closest=np.array(random.sample(range(MAX),2*N),dtype='uint32')
find_closest.sort()

filtered.tofile('/tmp/closest_index_after_test_filtered_forward.u32')
filtered.tofile('/tmp/closest_index_after_test_filtered_reverse.u32')
find_closest.tofile('/tmp/closest_index_after_test_find_closest_forward.u32')
find_closest.tofile('/tmp/closest_index_after_test_find_closest_reverse.u32')
