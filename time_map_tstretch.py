__doc__="""
This takes a sequence of time pairs describing a mapping between an input and an
output time effected by a time-stretching algorithm.
This provides functions that:
    - filter out impossible time mappings
    - give a list of input and output frames which can be used by
      time-stretching algorithms

General assumputions:
    H > 0
    L_w > 0

"""


class io_time_pair:
    def __init__(self,in_time,out_time):
        self.in_time=in_time
        self.out_time=out_time

def io_time_pair_in_time(x):
    return x.in_time

def io_time_pair_out_time(x):
    return x.out_time

def check_io_time_pairs_in_bounds(io_time_pairs,H,L_w,L_lb,L_ub):
    """
    Given a hop size H and window length L_w, determine if:
        the io_time_pair with the smallest in_time will map to an analysis frame
        with a start time >= L_lb (lower bound) AND
        the io_time_pair with the greatest out_time will map to analysis frame
        with start time <= L_ub - L_w (upper bound - frame length)
    """
    min_io_time_pair=sorted(io_time_pairs,key=io_time_pair_in_time)[0]
    max_io_time_pair=sorted(io_time_pairs,key=io_time_pair_in_time,reverse=True)[0]
    ret = True
    ret &= (min_io_time_pair.in_time - (min_io_time_pair.out_time % H)) >= L_lb
    ret &= (max_io_time_pair.in_time - (max_io_time_pair.out_time % H) + L) <= L_ub
    return ret

def filter_close_io_time_maps(
    # a list of [(intime,outtime),...] pairs
    io_time_pairs,
    # minimum allowed time between input times (typically the length of the analysis window L_w)
    t_in_min,
    # minimum allowed time between output times (typically H+L_w where H is the (output) hop size)
    t_out_min):
    """filter out impossible time mappings"""
    if len(io_time_pairs) < 1:
        return ValueError('io_time_pairs must contain at least 1 pair')
    assert(t_in_min>0)
    assert(t_out_min>0)
    #n_in=1
    n_out=1
    # in the case you wanted to preallocate this array, it would have a maximum
    # length the same as the length of io_time_pairs
    filtered_time_pairs=[io_time_pairs[0]]
    for n_in in range(1,len(io_time_pairs)):
        in_diff = io_time_pairs[n_in].in_time - filtered_time_pairs[n_out-1].in_time
        out_diff = io_time_pairs[n_in].out_time - filtered_time_pairs[n_out-1].out_time
        if (in_diff >= t_in_min) and (out_diff >= t_out_min):
            filtered_time_pairs.append(io_time_pairs[n_in])
            n_out+=1
    return filtered_time_pairs

class locked_frame_info_out:
    def __init__(self,
        # the start time of the frame
        start_time,
        # the pair contained in this frame, whose out time must be at its time
        # in the output
        # the output time is absolute (relative to the beginning of the time
        # series (recording))
        locked_time_pair,
        # the number of frames to the previous locked frame, if frames are
        # indexed sequentially, then the distance between 2 frames is
        # frame_index_1 - frame_index_0 where the convention is that a distance
        # is positive when frame_index_1 refers to a frame later than
        # frame_index_0.
        nframes_to_prev_locked):
        self.start_time=start_time
        self.locked_time_pair=locked_time_pair
        self.nframes_to_prev_locked=nframes_to_prev_locked
        
def out_locked_frames(
    # a sequence of filtered pairs describing an input and output mapping
    io_time_pairs,
    # the output (synthesis) hop size
    H,
    # the output (synthesis) window length
    L_w):
    # the first frame is always locked
    # in the case you wanted to preallocate "frames", it needs to have length of
    # at least len(io_time_pairs)*ceil(L_w/H)
    frames=[]
    # n_frame=0
    # sample number
    n=0
    # distance between frames in number of frames
    frame_dist=0
    for pair in io_time_pairs:
        while n <= pair.out_time:
            if (n+L_w) > pair.out_time:
                frames.append(n,pair,frame_dist)
                frame_dist = 0
                # n_frame += 1
            n += H
            frame_dist += 1
    return frames

def in_locked_frames(
    # the output frames (e.g., as determined by out_locked_frames)
    out_frames):
    """ Returns a list of sample numbers representing the start times of the
    input (analysis frames)"""
    assert(len(out_frames)>0)
    # in case you wanted to preallocate in_frames, it needs to have size
    # (sum([f.frame_dist for f in out_frames])+1)
    in_frames=[out_frames[0].locked_time_pair.in_time
                -(out_frames[0].locked_time_pair.out_time
                    - out_frames[0].start_time)]
    n_out_frames = len(out_frames)
    # n_out_frame = 1
    n_in_frame = 1
    for out_frame in out_frames:
        # the distance from the start of the in_frame to the in_time should be
        # the same as the distance between the start of the out_frame and the
        # out_time
        locked_frame_start_time = out_frame.locked_time_pair.in_time
                -(out_frame.locked_time_pair.out_time
                    - out_frame.start_time)
        local_hop_size = (locked_frame_start_time -
            in_frames[n_in_frame-1])/out_frame.nframes_to_prev_locked
        h_local = in_frames[n_in_frame-1]
        for n in range(out_frame.nframes_to_prev_locked-1):
            # in_frames[n_in_frame] = round(h_local)
            in_frames.append(round(h_local))
            h_local+=local_hop_size
            n_in_frame += 1
        in_frames.append(locked_frame_start_time)
        n_in_frame += 1
    return in_frames
