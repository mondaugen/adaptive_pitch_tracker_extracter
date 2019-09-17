from time_map_tstretch import *
from time_map_tstretch_test_common import *
from functools import reduce

L_w=8
tests_in_locked_frames=in_locked_frames(locked_out_frames_H4_W8_filtered_W8_time_pairs,L_w)
for i,ilf in enumerate(tests_in_locked_frames):
    print(i,ilf)
assert(reduce(lambda x,y: x and y,[x == y for x,y in
zip(tests_in_locked_frames,in_locked_frames_H4_W8_filtered_W8)]))
