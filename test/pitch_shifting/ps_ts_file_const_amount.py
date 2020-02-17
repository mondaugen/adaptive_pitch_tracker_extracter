# Simply pitch shift and or time-stretch a file by a constant amount.
# Detects the attacks and preserves them.

# Gate signal is just on for the full file and ADSR makes a very slight ramp at
# the beginning and ends of the signal.
# At the attacks, the time stretcher is reset and played for the time required
# to fill the envelope
# The time-stretcher uses the audio-rate time-stretch and pitch-shift signals
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
import lfo

# Environment variables
SR=common.get_env('SR',default=16000,conv=int)
W=common.get_env('W',default=1024,conv=int)
H=common.get_env('H',default=256,conv=int)
IN_FILE=common.get_env('IN_FILE',default='/tmp/in.f64')
OUT_FILE=common.get_env('OUT_FILE',default='/tmp/out.f64')
LMAX_FILT_RATE=common.get_env('LMAX_FILT_RATE',default=SR,conv=float)
# Gate parameters
ADSR_ATTACK=common.get_env('ADSR_ATTACK',default=0.1,conv=float)
ADSR_RELEASE=common.get_env('ADSR_RELEASE',default=0.1,conv=float)
# Time-stretch LFO parameters
TS=common.get_env('TS',default=1,conv=float)
# Pitch-shift LFO parameters
PS=common.get_env('PS',default=1,conv=float)

# get signal and adjust length
x=np.fromfile(IN_FILE,dtype='float64')
N=0
while N < len(x):
    N+=H
x=np.concatenate((x,np.zeros(N-len(x))))
# add dither
x+=np.random.standard_normal(len(x))*1e-8

# make time-stretch signal, just constant for whole file
ts_sig=np.ones(N)*TS
# make pitch-shift signal, just constant for whole file
ps_sig=np.ones(N)*PS
# make position signal, just 0 for whole file
pos_sig=np.zeros(N)

# estimate attack times
attack_time_pairs=attack_finder.attacks_from_spectral_diff(x,lmax_filt_rate=LMAX_FILT_RATE)
attack_times=np.array([b for a,b in attack_time_pairs])
# make attack avoider
#av=time_map_tstretch.attack_avoider(attack_times,-H,H+W,H)
av=time_map_tstretch.attack_avoider(
    attack_times,
    -3*H,
    W+2*H,
    H)

# synthesize gate signal, just on the whole time
gate_sig=np.ones(N)#np.concatenate((np.zeros(1),np.ones(N-2),np.zeros(1)))

rs=envelopes.region_segmenter(H)
# make way to look up signal (the following is fast at the end-points)
wl=window_tools.windowed_lookup(x,W)
pv=pvoc_synth(
    signal.get_window('hann',W),
    signal.get_window('hann',W),
    W,
    H,
    lambda n: wl.access(n))
aaa=time_map_tstretch.attack_avoid_access(
    lambda t,r: pv.process(int(np.round(t)),r),
    av)
ps=pitch_shift.pitch_shifter(
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
                ps.set_pos_at_block_start(pos_sig[n_+s])
            y[n_+s:n_+e]=ps.process(ps_sig[n_+s:n_+e],ts_sig[n_+s:n_+e])*en['adsr'][s:e]

y.tofile(OUT_FILE)
