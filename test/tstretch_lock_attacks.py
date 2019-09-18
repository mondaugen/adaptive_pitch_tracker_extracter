__doc__="""
take a soundfile and the score used to render it (see scripts/synth_af_score.py)
and use this score to determine the points that are fixed in a time stretch.
Produce a file that is the time-stretched version of the input file, but with
the attacks preserved.
"""

import classic_puckette_timestretch
import synth_af_score
import common
import time_map_tstretch
import numpy as np

if __name__ == '__main__':
    SCORE_FILE=common.get_env('SCORE_FILE',check_if_none=True)
    SOUND_FILE=common.get_env('SOUND_FILE',check_if_none=True)
    OUTPUT=common.get_env('OUTPUT',check_if_none=True)
    # smaller numbers mean longer file
    TSTRETCH_FACTOR=common.get_env('TSTRETCH_FACTOR',default=1,conv=float)
    # move the locked points forward or backward in time by a fixed amount
    LOCK_OFFSET=common.get_env('LOCK_OFFSET',default=0,conv=int)
    # a time after the lock time within which all the points are locked
    # this value is multiplied by H+W because we cannot lock 2 points closer
    # than that amount of time
    LOCK_WINDOW=common.get_env('LOCK_WINDOW',default=0,conv=int)

    H=common.get_env('H',default=256,conv=int)
    W=common.get_env('W',default=1024,conv=int)
    SAMPLE_RATE=common.get_env('SAMPLE_RATE',16000,float)
    x=np.fromfile(SOUND_FILE,'float64')
    orig_len_x=len(x)
    # instead of adding LOCK_OFFSET to the input points, we prepend zeros to the
    # input file so that reads before the beginning of the file can happen. This
    # will also offset the locked points.
    locked_input_points=synth_af_score.extract_ts_samps(SCORE_FILE,SAMPLE_RATE)
    # add point near end because we want to play to the end
    locked_input_points+=[orig_len_x-W]
    if LOCK_OFFSET > 0:
        x=np.concatenate((x,np.random.standard_normal(LOCK_OFFSET).astype(x.dtype)*1e-8))
    elif LOCK_OFFSET < 0:
        x=np.concatenate((np.random.standard_normal(-LOCK_OFFSET).astype(x.dtype)*1e-8,x))
        
    io_time_pairs=[time_map_tstretch.io_time_pair(ip,ip/TSTRETCH_FACTOR) for ip
        in locked_input_points]
    if LOCK_WINDOW > 0:
        augmented_io_time_pairs=[]
        lw=(H+W)*LOCK_WINDOW
        for pair in io_time_pairs:
            augmented_io_time_pairs.append(pair)
            if (pair.in_time + lw) < (orig_len_x-W):
                # we don't want to add any points that would make an analysis
                # window extend beyond the end of the sound file
                augmented_io_time_pairs.append(
                    time_map_tstretch.io_time_pair(
                        pair.in_time+lw,
                        pair.out_time+lw))
        io_time_pairs=augmented_io_time_pairs
    # sort because the lock window could have introduced points extending past the next points
    io_time_pairs=time_map_tstretch.sort_io_time_pairs_in_time(io_time_pairs)
    assert(time_map_tstretch.check_io_time_pairs_in_bounds(io_time_pairs,H,W,0,len(x)))
    # we also want to map 0 -> 0. This will be okay (it will not cause an
    # out-of-bounds reading) because the time stretcher doesn't use a previous
    # input frame to synthesize the first output frame
    #io_time_pairs=[time_map_tstretch.io_time_pair(0,0)]+io_time_pairs
    io_time_pairs=time_map_tstretch.filter_close_inc_io_time_maps(io_time_pairs,H+W,H+W)
    out_locked_frames=time_map_tstretch.out_locked_frames(io_time_pairs,H,W)
    in_locked_frames=time_map_tstretch.in_locked_frames(out_locked_frames,W)
    in_locked_frames=np.array(in_locked_frames,dtype='int')
    print(in_locked_frames)
    y=classic_puckette_timestretch.time_stretch_arb_times(x,in_locked_frames,H,W)
    y=common.normalize(y)
    y.tofile(OUTPUT)

    
