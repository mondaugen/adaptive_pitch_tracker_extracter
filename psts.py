# Common pitch-shifting and time-stretching routines

import numpy as np
from scipy import signal
import pitch_shift
import window_tools
from classic_puckette_timestretch import pvoc_synth
import matplotlib.pyplot as plt
import rel_del_line
import time_map_tstretch
import attack_finder
import common
import envelopes
import math

def adjust_x_for_time_stretch(x,TS,H):
    new_len_=math.ceil(len(x)/TS)
    new_len = 0
    while new_len < new_len_:
        new_len += H
    if new_len < len(x):
        return (x[:new_len],new_len)
    elif new_len > len(x):
        len_diff = new_len - len(x)
        return (np.concatenate((x,np.zeros(len_diff,dtype=x.dtype))),new_len)
    return (x,new_len)

def pad_x_and_dither(x,H,amt=1e-6):
    N=0
    while N < len(x):
        N+=H
    x=np.concatenate((x,np.zeros(N-len(x))))
    # add dither
    x+=np.random.standard_normal(len(x))*amt
    return (x,N)

# you need to filter out attack times that are closer than awin_len + 1
def psts_const_amount(
    x, # input signal
    # attack times, e.g., estimate with 
    # attack_time_pairs=attack_finder.attacks_from_spectral_diff(x,
    # lmax_filt_rate=LMAX_FILT_RATE)
    # attack_times=np.array([b for a,b in attack_time_pairs])
    attack_times,
    # see time_map_tstretch.attack_avoider for more information about awin_start
    # and awin_len
    awin_start=None,
    awin_len=None,
    # type of analysis window
    awin_type='hann',
    # type of synthesis window
    swin_type='hann',
    SR=16000,
    W=1024,
    H=256,
    # Gate parameters
    ADSR_ATTACK=0.1,
    ADSR_RELEASE=0.1,
    # Time-stretch LFO parameters
    TS=1,
    # Pitch-shift LFO parameters
    PS=1,
    # record the look-up times in this file
    look_up_times_path='/tmp/look_up_times.u32',
    # record the reset times in this file
    reset_times_path='/tmp/reset_times.u32',
    ):

    if awin_start is None:
        awin_start=-3*H
    if awin_len is None:
        awin_len=W+2*H

    N=0
    while N < len(x):
        N+=H
    x=np.concatenate((x,np.zeros(N-len(x))))
    x,N=adjust_x_for_time_stretch(x,TS,H)
    # add dither
    x+=np.random.standard_normal(len(x))*1e-8

    # make time-stretch signal, just constant for whole file
    ts_sig=np.ones(N)*TS
    # make pitch-shift signal, just constant for whole file
    ps_sig=np.ones(N)*PS
    # make position signal, just 0 for whole file
    pos_sig=np.zeros(N)

    # make attack avoider
    av=time_map_tstretch.attack_avoider(
        attack_times,
        awin_start,
        awin_len,
        H)

    # synthesize gate signal, just on the whole time
    gate_sig=np.ones(N)

    rs=envelopes.region_segmenter(H)
    # make way to look up signal (the following is fast at the end-points)
    wl=window_tools.windowed_lookup(x,W)
    pv=pvoc_synth(
        signal.get_window(awin_type,W),
        signal.get_window(swin_type,W),
        W,
        H,
        lambda n: wl.access(n))

    # make way to get samples at given time from pvoc
    look_up_times_fd=open(look_up_times_path,'w')
    reset_times_fd=open(reset_times_path,'w')
    def _pvoc_get_samples(t,r):
        t=int(np.round(t))
        # record the look up times
        np.array([t],dtype='uint32').tofile(look_up_times_fd)
        # if reset occurred, record this too
        if r:
            np.array([t],dtype='uint32').tofile(reset_times_fd)
        return pv.process(t,r)

    aaa=time_map_tstretch.attack_avoid_access(
        _pvoc_get_samples,
        av)
    psps=pitch_shift.pitch_shifter(
    aaa,
    B=H)
    adsr=envelopes.gate_to_adsr(ADSR_ATTACK*SR,1,1,ADSR_RELEASE*SR)

    y=np.zeros_like(x)
    for n_ in range(0,N,H):
        H_=min(N-n_,H)
        en,st,ed,ac,sta=adsr.gate_to_adsr_env_start_end_active(gate_sig[n_:n_+H_])
        ant,ans,ane,regs=rs.region_segmenter_update(st,ed,ac)
        for s,e in regs:
            if e>s:
                l=e-s
                if st[s] > 0:
                    pv.reset_past_output()
                    psps.set_pos_at_block_start(pos_sig[n_+s])
                y[n_+s:n_+e]=psps.process(
                ps_sig[n_+s:n_+e],ts_sig[n_+s:n_+e])*en['adsr'][s:e]

    look_up_times_fd.close()
    reset_times_fd.close()

    return y

# pitch shift in real-time, varying the pitch shift
def ps_lfo_rt(
    x, # input signal
    ps, # pitch shift amounts
    # time in samples for max filtering IR to fall to 0.01 of original value
    lmax_filt_rate=16000,
    # minimum threshold of local power for local maximum in attack estimation
    # signal (spectral difference) to be accepted
    ng_th=-60,
    # type of analysis window
    awin_type='hann',
    # type of synthesis window
    swin_type='hann',
    SR=16000,
    W=1024,
    H=256,
    # Gate parameters
    ADSR_ATTACK=0.1,
    ADSR_RELEASE=0.1,
    # for debugging and A/B comparison, if true, no compensation is made for
    # attacks
    always_ignore_attack=False):

    N=len(x)
    if (len(ps) != N) or ((N % H) != 0):
        raise ValueError

    # make time-stretch signal, just constant for whole file
    ts_sig=np.ones(N)
    # make position signal, just 0 for whole file
    pos_sig=np.zeros(N)

    # convert lmax_filt_rate in samples to hops
    lmax_filt_rate_h=int(np.ceil(lmax_filt_rate/H))

    # the size of the safe region in which an analysis window can sit without
    # covering any attacks
    #R=2*M+W+H
    # attacks can come no closer than this many hops
    attack_freq_limit=1#int(np.ceil((R+1)/H))

    # analysis object for attacks
    afsd=attack_finder.attacks_from_spectral_diff_rt(
    W=W,
    H=H,
    lmax_filt_rate=lmax_filt_rate_h,
    attack_freq_limit=attack_freq_limit,
    ng_th=ng_th)

    # synthesize gate signal, just on the whole time
    gate_sig=np.ones(N)

    rs=envelopes.region_segmenter(H)

    Pmin=np.min(ps)
    Pmax=np.max(ps)
    # number of extra values requested by interpolator
    I=pitch_shift.default_get_interpolator_n_points(0)
    # max number of reads per write when getting values for pitch shifting
    Nr=int(np.ceil((Pmax*H+I)/H))
    # max number of writes per read when getting values for pitch shifting
    Nw=int(np.ceil(1/Pmin))

    rtpaac=time_map_tstretch.real_time_ps_attack_avoid_controller(W,H,Nr,Nw)
    def _rtaac_read():
        samps,reset=rtpaac.read()
        return (samps,reset)
    pv=pvoc_synth(
        signal.get_window(awin_type,W),
        signal.get_window(swin_type,W),
        W,
        H,
        get_samples_no_time=_rtaac_read)
    psps=pitch_shift.pitch_shifter(
    # Actually time and reset are ignored because we've specified
    # get_samples_no_time
    lambda t: pv.process(t,False),
    B=H)
    adsr=envelopes.gate_to_adsr(ADSR_ATTACK*SR,1,1,ADSR_RELEASE*SR)

    y=np.zeros_like(x)
    for n_ in range(0,N,H):
        H_=min(N-n_,H)
        en,st,ed,ac,sta=adsr.gate_to_adsr_env_start_end_active(gate_sig[n_:n_+H_])
        ant,ans,ane,regs=rs.region_segmenter_update(st,ed,ac)
        for s,e in regs:
            if e>s:
                l=e-s
                if st[s] > 0:
                    pv.reset_past_output()
                    psps.set_pos_at_block_start(pos_sig[n_+s])
                samps=x[n_+s:n_+e]
                attack=afsd(samps)
                rtpaac.write(samps,(attack > 0)and(not always_ignore_attack))
                y[n_+s:n_+e]=psps.process(
                ps[n_+s:n_+e],ts_sig[n_+s:n_+e])*en['adsr'][s:e]
    return y


# uses real-time method of attack estimation and avoidance
def psts_const_amount_rt(
    x, # input signal
    # time in samples for max filtering IR to fall to 0.01 of original value
    lmax_filt_rate=16000,
    # minimum threshold of local power for local maximum in attack estimation
    # signal (spectral difference) to be accepted
    ng_th=-60,
    # type of analysis window
    awin_type='hann',
    # type of synthesis window
    swin_type='hann',
    SR=16000,
    W=1024,
    H=256,
    # Gate parameters
    ADSR_ATTACK=0.1,
    ADSR_RELEASE=0.1,
    # Time-stretch not possible in real-time
    # Pitch-shift LFO parameters
    PS=1,
    # for debugging and A/B comparison, if true, no compensation is made for
    # attacks
    always_ignore_attack=False):

    x,N=pad_x_and_dither(x,H,1e-8)

    # make pitch-shift signal, just constant for whole file
    ps=np.ones(N)*PS

    return ps_lfo_rt(
        x,
        ps,
        lmax_filt_rate=lmax_filt_rate,
        ng_th=ng_th,
        awin_type=awin_type,
        swin_type=swin_type,
        SR=SR,
        W=W,
        H=H,
        ADSR_ATTACK=ADSR_ATTACK,
        ADSR_RELEASE=ADSR_RELEASE,
        always_ignore_attack=always_ignore_attack)
