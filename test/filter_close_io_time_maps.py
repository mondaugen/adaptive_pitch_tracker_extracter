from time_map_tstretch import *
from time_map_tstretch_test_common import *
from functools import reduce

L_w=8
tests_filtered_maps=filter_close_inc_io_time_maps(io_time_pairs,L_w,L_w)
assert(reduce(lambda x,y: x and y,[x == y for x,y in zip(tests_filtered_maps,io_filtered_W8_time_pairs)]))
