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

from functools import reduce
import numpy as np

class io_time_pair:
    def __init__(self,in_time,out_time):
        self.in_time=in_time
        self.out_time=out_time
    def __eq__(self,other):
        return (self.in_time==other.in_time) and (self.out_time==other.out_time)
    def __str__(self):
        return '(i: %d, o: %d)' % (self.in_time,self.out_time)

def io_time_pair_in_time(x):
    return x.in_time

def sort_io_time_pairs_in_time(io_time_pairs,reverse=False):
    return sorted(io_time_pairs,key=io_time_pair_in_time,reverse=reverse)

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
    min_io_time_pair=sort_io_time_pairs_in_time(io_time_pairs)[0]
    max_io_time_pair=sort_io_time_pairs_in_time(io_time_pairs,reverse=True)[0]
    ret = True
    # for the first output frame, the locked feature is mapped onto the first
    # sample of the output frame, so the locked feature is also on the first
    # sample of the input frame, so nothing needs to be added to the in_time to
    # get the frame start time
    ret &= min_io_time_pair.in_time >= L_lb
    # we will set the first output frame's start time to be equal to the output
    # time of the locked feature. This means there are
    # (max_io_time_pair.out_time - min_io_time_pair.out_time) / H frames (hops)
    # until the last last locked feature.  Because we want the locked output
    # feature to map to the same sample number of the input frame (i.e., we want
    # it to have the same distance from the input frame's start as it has from
    # the output frame's start) we subtract this remainder multiplied by the hop
    # size from min_io_time_pair.in_time to get the last input frame's start
    # time. The length of the input frame is L, so the last valid sample is
    # input_frame.start_time + L - 1 (indexing starts at 0), so we check that
    # this value <= the upper bound
    ret &= (max_io_time_pair.in_time
        - ((max_io_time_pair.out_time 
            - min_io_time_pair.out_time) % H) + L_w) <= L_ub
    return ret

def check_io_time_pairs_out_time_increasing(io_time_pairs):
    """ Check to see all the out_times are monotonically increasing, i.e.,
    io_time_pairs[i+1]>=io_time_pairs[i] for all i """
    return reduce(lambda x,y: x and y,
        [(x.out_time-y.out_time)>=0 for x,y in zip(io_time_pairs[1:],io_time_pairs[:-1])])

def filter_close_inc_io_time_maps(
    # a list of [(intime,outtime),...] pairs
    io_time_pairs,
    # minimum allowed time between input times (typically the length of the analysis window L_w)
    t_in_min,
    # minimum allowed time between output times (typically H+L_w where H is the (output) hop size)
    t_out_min):
    """
    filter out time mappings smaller than some values
    this function is designed to work on monotonically increasing mappings
    """
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
    def __eq__(self,other):
        return ((self.start_time==other.start_time)
            and (self.locked_time_pair == other.locked_time_pair)
            and (self.nframes_to_prev_locked == other.nframes_to_prev_locked))
    def __str__(self):
        return (('(s: %d, ' % self.start_time)
            + str(self.locked_time_pair)
            + (', n: %d)' % self.nframes_to_prev_locked))
        
def out_locked_frames(
    # a sequence of filtered pairs describing an input and output mapping
    io_time_pairs,
    # the output (synthesis) hop size
    H,
    # the output (synthesis) window length
    L_w):
    assert(len(io_time_pairs)>0)
    assert(H>0)
    assert(L_w>0)
    # the first frame is always locked
    # in the case you wanted to preallocate "frames", it needs to have length of
    # at least len(io_time_pairs)*ceil(L_w/H)
    # the first output frame is simply set so that the first sample in the frame
    # is the first output time we want to map to
    # the number of frames to the previous frame is 0 because it is the first frame
    frames=[locked_frame_info_out(io_time_pairs[0].out_time,io_time_pairs[0],0)]
    # n_frame=0
    # sample number
    n=io_time_pairs[0].out_time+H
    # distance between frames in number of frames
    frame_dist=1
    for pair in io_time_pairs[1:]:
        while n <= pair.out_time:
            if (n+L_w) > pair.out_time:
                frames.append(locked_frame_info_out(n,pair,frame_dist))
                frame_dist = 0
                # n_frame += 1
            n += H
            frame_dist += 1
    return frames

def in_locked_frames(
    # the output frames (e.g., as determined by out_locked_frames)
    out_frames,
    # the analysis frame length
    L_w):
    """ Returns a list of sample numbers representing the start times of the
    input frames (analysis frames)"""
    assert(len(out_frames)>0)
    # in case you wanted to preallocate in_frames, it needs to have size
    # (sum([f.frame_dist for f in out_frames[1:]])+1)
    in_frames=[out_frames[0].locked_time_pair.in_time
                -(out_frames[0].locked_time_pair.out_time
                    - out_frames[0].start_time)]
    n_out_frames = len(out_frames)
    # n_out_frame = 1
    n_in_frame = 1
    for out_frame in out_frames[1:]:
        # the distance from the start of the in_frame to the in_time should be
        # the same as the distance between the start of the out_frame and the
        # out_time
        lock_input_frame_start_time = (out_frame.locked_time_pair.in_time
                -(out_frame.locked_time_pair.out_time
                    - out_frame.start_time))
        # the previous locked frame input frame should not contain the locked
        # point, and is therefore set to have a start time equal to the in_time
        # - the analysis frame length (so that the sample just after the last
        # sample in the previous locked frame is the locked point)
        # If the input points are filtered so that no input points are closer
        # than a window length, this input frame will always come at a time =>
        # the previous input frame (in_frames[n_in_frame-1])
        # otherwise, the local_hop_size could be negative, in which case the
        # in_frames might advance backwards, which is not really wrong, just
        # maybe surprising
        prev_lock_input_frame_start_time = out_frame.locked_time_pair.in_time - L_w
        if out_frame.nframes_to_prev_locked > 1:
            local_hop_size = (prev_lock_input_frame_start_time -
                in_frames[n_in_frame-1])/(out_frame.nframes_to_prev_locked-1)
            h_local = in_frames[n_in_frame-1]+local_hop_size
            for n in range(out_frame.nframes_to_prev_locked-1):
                # in_frames[n_in_frame] = round(h_local)
                in_frames.append(round(h_local))
                h_local+=local_hop_size
                n_in_frame += 1
        in_frames.append(lock_input_frame_start_time)
        n_in_frame += 1
    return in_frames

def get_contained_attack_times(attack_times, start, end):
    return attack_times[np.where((attack_times >= start)&(attack_times<end))[0]]

def contains_atime(attack_times, start, end):
    return len(get_contained_attack_times(attack_times,start,end)) > 0

class attack_avoider:
    def __init__(self,
        # Some data structure containing points in time that we want to preserve. It
        # shouldn't contain time points closer than awin_len, otherwise this
        # algorithm won't be able put a window in a place not containing one of
        # those time points.
        attack_times,
        # This is the time we want to put the analysis window. A possible/recommended
        # value for this is a number such that analysis_time + awin_start gives
        # the start time of the window, but other values might end up sounding better.
        awin_start,
        # A number such that analysis_time + awin_start + awin_len gives the first
        # sample just after the end of a window.
        awin_len,
        # a synthesis hop size
        H):
        self.attack_times = np.array(attack_times)
        self.awin_start = awin_start
        self.awin_len = awin_len
        self.H = H
        self.last_contained_atime = -1
        self.last_output_time = -1
        self.safe = True
    def adjust(self, analysis_time):
        """ Adjust this analysis_time so it doesn't land on an attack.  Returns
        a tuple (adjusted_time,reset) where adjusted_time is the time adjusted
        if it would smear an attack and reset is a boolean which is true if the
        returned time is the first frame containing an attack. This boolean can
        be used to know when a time stretching algorithm should reset the phase
        of the output (to avoid catastrophic smearing). """
        contained_atime = get_contained_attack_times(
            self.attack_times,
            analysis_time + self.awin_start,
            analysis_time + self.awin_start + self.awin_len)
        # initially we propose the analysis_time but check to see if it would
        # smear the preserved attack time
        ret = analysis_time
        reset=False
        if len(contained_atime) > 0:
            contained_atime = contained_atime[0]
            if contained_atime == self.last_contained_atime:
                # If contained_atime is the same as last time, output a window a
                # normal hop size's distance from the last output (if the
                # analysis window isn't already safe)
                ret = self.last_output_time
                if not self.safe:
                    ret += self.H
            else:
                # use the new contained atime and store it
                self.last_contained_atime = contained_atime
                ret = analysis_time
                reset=True
        self.last_output_time = ret
        if contains_atime(
            self.attack_times,
            self.last_output_time + self.awin_start,
            self.last_output_time + self.awin_start + self.awin_len):
            # we are not in a safe position
            self.safe = False
        else:
            self.safe = True
        return (ret,reset)

class attack_avoid_access:
    """
    Many routines need a way to get samples at a particular time. This callable
    object (like a function) provides n samples requested at a time t. This
    adjusts the requested times to keep attacks from getting "smeared" if an
    attack_avoider is registered with it.
    """
    def __init__(self,
        # a function accepting a time t and a flag r that when true indicates
        # that a reset should be performed. r will be true when the time passed
        # is the first of a series of adjusted times
        get_samples,
        av):
        self.av=av
        self.get_samples=get_samples
    def __call__(self,t):
        atime,reset=self.av.adjust(int(np.round(t)))
        return self.get_samples(atime,reset)
