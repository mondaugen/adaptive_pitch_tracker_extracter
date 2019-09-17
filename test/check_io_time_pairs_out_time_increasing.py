from time_map_tstretch import *
from time_map_tstretch_test_common import *

bad_io_time_pairs=io_time_pairs[len(io_time_pairs)//2:]+io_time_pairs[:len(io_time_pairs)//2]
assert(check_io_time_pairs_out_time_increasing(io_time_pairs)==True)
assert(check_io_time_pairs_out_time_increasing(bad_io_time_pairs)==False)
