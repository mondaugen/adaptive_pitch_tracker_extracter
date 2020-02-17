# A gate signal is converted to an ADSR signal
# At the attacks, the time stretcher is reset and played for the time required
# to fill the envelope
# The time-stretcher uses the audio-rate time-stretch and pitch-shift signals
import numpy as np
from scipy import signal
from pitch_shift import pitch_shifter
import window_tools
from classic_puckette_timestretch import pvoc_synth
import matplotlib.pyplot as plt
import rel_del_line
from time_map_tstretch import attack_avoider
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
GATE_F0=common.get_env('GATE_F0',default=1,conv=float)/SR
GATE_F1=common.get_env('GATE_F1',default=GATE_F0*SR,conv=float)/SR
# Time-stretch LFO parameters
TS_FM0=common.get_env('TS_FM0',default=0.5,conv=float)/SR
TS_FM1=common.get_env('TS_FM1',default=TS_FM0*SR,conv=float)/SR
TS_MIN=common.get_env('TS_MIN',default=0.5,conv=float)
TS_MAX=common.get_env('TS_MAX',default=1.5,conv=float)
# Pitch-shift LFO parameters
PS_FM0=common.get_env('PS_FM0',default=0.75,conv=float)/SR
PS_FM1=common.get_env('PS_FM1',default=PS_FM0*SR,conv=float)/SR
PS_MIN=common.get_env('PS_MIN',default=0.5,conv=float)
PS_MAX=common.get_env('PS_MAX',default=1.5,conv=float)
# Position signal LFO parameters
POS_FM0=common.get_env('POS_FM0',default=0.75,conv=float)/SR
POS_FM1=common.get_env('POS_FM1',default=POS_FM0*SR,conv=float)/SR
POS_MIN=max(common.get_env('POS_MIN',default=0,conv=float),0)
POS_MAX=min(common.get_env('POS_MAX',default=1,conv=float),1)

# get signal and adjust length
x=np.fromfile(IN_FILE,dtype='float64')
N=0
while N < len(x):
    N+=H
x=np.concatenate((x,np.zeros(N-len(x))))
# add dither
x+=np.random.standard_normal(len(x))*1e-8

# make time-stretch signal
ts_sig=lfo.chirp(N,TS_FM0,TS_FM1,x_min=TS_MIN,x_max=TS_MAX).sinusoid()
# make pitch-shift signal
ps_sig=lfo.chirp(N,PS_FM0,PS_FM1,x_min=PS_MIN,x_max=PS_MAX).sinusoid()
# make position signal
pos_sig=lfo.chirp(N,POS_FM0,POS_FM1,x_min=POS_MIN,x_max=POS_MAX).sinusoid()

# make way to look up signal (the following is fast at the end-points)
wl=window_tools.windowed_lookup(x,W)
def wl_access(t,l):
    return wl.access(t)

# make the phase vocoder synthesizer, which uses the wl_access function to look
# up a window's worth of samples
pv=pvoc_synth(
    signal.get_window('hann',W),
    signal.get_window('hann',W),
    W,
    H,
    wl_access)

# estimate attack times
attack_time_pairs=attack_finder.attacks_from_spectral_diff(x,lmax_filt_rate=LMAX_FILT_RATE)
attack_times=np.array([b for a,b in attack_time_pairs])
# make attack avoider
av=attack_avoider(attack_times,-H,H+W,H)

# synthesize gate signal
gate_sig=lfo.chirp(N,GATE_F0,GATE_F1).squarewave()

# make ADSR object
adsr=envelopes.gate_to_adsr(SR*0.1, SR*0.1, 0.5, SR*.1)

# make region segmenter
rs=envelopes.region_segmenter(H)

class time_stretch_ts_access:
    """
    The pitch shifter needs a way to get samples to fill a ring buffer, which it
    then resamples to effect the pitch shift. This callable object (like a
    function) provides n samples requested at a time t. This particular one uses a
    phase vocoder to do it in a way so that the phase of the signal is
    smooth from one block to the next, even if the two consecutive times
    looked up are not separated by a hop size.
    This also adjusts the requested times to keep attacks from getting
    "smeared" if an attack_avoider is registered with it.
    """
    def __init__(self,init_s,pv,av=None):
        self.pv=pv
        self.s=init_s
        self.t_off=0
        self.av=av
    def __call__(self,t,n):
        # this doesn't correctly time stretch the signal.
        # TODO: The goal is: for a given time-stretch and pitch-shift and start
        # position, we want the first sample to correspond to the same start
        # position were there a time-stretch and pitch-shift factor of 1
        t*=self.s
        t+=self.t_off
        if self.av is not None:
            atime,reset=self.av.adjust(int(np.round(t)))
        else:
            atime=int(np.round(t))
            reset=False
        return pv.process(atime,reset)

# make an instance of the time_stretch_ts_access class, to give to pitch_shifter
ts_access=time_stretch_ts_access(1,pv,av)
ps=pitch_shifter(ts_access,B=H)

y=np.zeros_like(x)
for n in np.arange(0,N-H,H):
    (adsr_dict,
    start,
    end,
    state,
    adsr_states)=adsr.gate_to_adsr_env_start_end_active(gate_sig[n:n+H])
    ant,ans,ane,regs=rs.region_segmenter_update(start,end,state)
    for s,e in regs:
        if e>s:
            l=e-s
            # set the time stretch amount by taking the mean of the stretch signal
            ts_access.s=np.mean(ts_sig[n+s:n+e])
            # hmm, not sure if this is correct...
            if start[s] > 0:
                ts_access.t_off=-ps.time_at_block_start+pos_sig[n+e]*len(x)
            y[n+s:n+e]=ps.process(ps_sig[n+s:n+e])*adsr_dict['adsr'][s:e]

y.tofile(OUT_FILE)
