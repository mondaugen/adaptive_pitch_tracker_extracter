from time_map_tstretch import *
from time_map_tstretch_test_common import *
from functools import reduce

L_w=8
H=4
tests_out_locked_frames=out_locked_frames(io_filtered_W8_time_pairs,H,L_w)
for x,y in zip(tests_out_locked_frames,
    locked_out_frames_H4_W8_filtered_W8_time_pairs):
    print("test: %s; true: %s; eq? %s\n" % (str(x),str(y),str(x==y)))
assert(reduce(lambda x,y: x and y,[x == y for x,y in zip(tests_out_locked_frames,locked_out_frames_H4_W8_filtered_W8_time_pairs)]))
