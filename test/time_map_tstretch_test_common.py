from time_map_tstretch import *

time_pairs=[
(7,0),
(11,8),
(19,25),
(30,30),
(39,48),
(47,56),
(64,64),
(80,72)
]
io_time_pairs=[io_time_pair(in_time,out_time) for in_time,out_time in time_pairs]

time_pairs_2=[
(7,1),
(11,8),
(19,25),
(30,30),
(39,48),
(47,56),
(64,64),
(80,72)
]
io_time_pairs_2=[io_time_pair(in_time,out_time) for in_time,out_time in time_pairs_2]

filtered_W8_time_pairs=[
(7,0),
(19,25),
(39,48),
(47,56),
(64,64),
(80,72)
]
io_filtered_W8_time_pairs=[io_time_pair(in_time,out_time) for in_time,out_time in filtered_W8_time_pairs]

out_frames_H4_W8_filtered_W8_time_pairs=[
(0,io_filtered_W8_time_pairs[0],0),
(20,io_filtered_W8_time_pairs[1],5),
(24,io_filtered_W8_time_pairs[1],1),
(44,io_filtered_W8_time_pairs[2],5),
(48,io_filtered_W8_time_pairs[2],1),
(52,io_filtered_W8_time_pairs[3],1),
(56,io_filtered_W8_time_pairs[3],1),
(60,io_filtered_W8_time_pairs[4],1),
(64,io_filtered_W8_time_pairs[4],1),
(68,io_filtered_W8_time_pairs[5],1),
(72,io_filtered_W8_time_pairs[5],1)]
locked_out_frames_H4_W8_filtered_W8_time_pairs=[locked_frame_info_out(*tup) for
    tup in out_frames_H4_W8_filtered_W8_time_pairs]

in_locked_frames_H4_W8_filtered_W8=[
7,
8,
9,
10,
11,
14,
18,
21,
24,
28,
31,
35,
39,
43,
47,
60,
64,
76,
80]
