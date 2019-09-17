from time_map_tstretch import *
from time_map_tstretch_test_common import *

assert(check_io_time_pairs_in_bounds(io_time_pairs,4,8,0,88)==True)
assert(check_io_time_pairs_in_bounds(io_time_pairs,4,8,0,87)==False)
assert(check_io_time_pairs_in_bounds(io_time_pairs_2,4,8,0,87)==True)
