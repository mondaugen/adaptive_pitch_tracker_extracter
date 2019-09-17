from time_map_tstretch import *
from time_map_tstretch_test_common import *
from functools import reduce

tests_in_locked_frames=in_locked_frames(locked_out_frames_H4_W8_filtered_W8_time_pairs)
for i,ilf in enumerate(tests_in_locked_frames):
    print(i,ilf)
