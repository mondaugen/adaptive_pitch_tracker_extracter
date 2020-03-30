# tests real-time attack avoid and pitch shift

import psts
import common
import numpy as np
import attack_finder
import lfo

IN_FILE=common.get_env('IN_FILE',default='/tmp/in.f64')
OUT_FILE=common.get_env('OUT_FILE',default='/tmp/out.f64')
SR=common.get_env('SR',default=16000,conv=int)
LMAX_FILT_RATE=common.get_env('LMAX_FILT_RATE',default=SR,conv=float)
N_CHANS=common.get_env('N_CHANS',default=1,conv=int)
W=common.get_env('W',default=1024,conv=int)
H=common.get_env('H',default=256,conv=int)
M=common.get_env('M',default=W,conv=int)
MIN_ATTACK_DIST=common.get_env('MIN_ATTACK_DIST',default=W+2*H,conv=int)
ATTACK_EST_CHAN=common.get_env('ATTACK_EST_CHAN',default=0,conv=int)
ATTACK_EST_W=common.get_env('ATTACK_EST_W',default=W,conv=int)
ATTACK_EST_H=common.get_env('ATTACK_EST_H',default=H,conv=int)
NG_TH=common.get_env('NG_TH',default=-60,conv=float)
ALWAYS_IGNORE_ATTACKS=common.get_env('ALWAYS_IGNORE_ATTACKS',default=0,conv=int)
F0=common.get_env('F0',default=1,conv=float)
F1=common.get_env('F1',default=1,conv=float)
PMIN=common.get_env('PMIN',default=1,conv=float)
PMAX=common.get_env('PMAX',default=1,conv=float)
PHASE=common.get_env('PHASE',default=0,conv=float)
LFO_TYPE=common.get_env('LFO_TYPE',default='sinusoid')

realtime=True

x=np.fromfile(IN_FILE,dtype='float64')
x,N=psts.pad_x_and_dither(x,N_CHANS*H)
N/=N_CHANS
ps_maker=lfo.chirp(
        N,
        F0/SR,
        F1/SR,
        x_min=PMIN,
        x_max=PMAX,
        phase=PHASE)

y_sigs=[]

for chan in range(N_CHANS):
    y_sigs.append(psts.ps_lfo_rt(
        x[chan::N_CHANS],
        eval("ps_maker.%s()" % (LFO_TYPE,)),
        lmax_filt_rate=LMAX_FILT_RATE,
        SR=SR,
        W=W,
        H=H,
        ADSR_ATTACK=common.get_env('ADSR_ATTACK',default=0.1,conv=float),
        ADSR_RELEASE=common.get_env('ADSR_RELEASE',default=0.1,conv=float),
        always_ignore_attack=bool(ALWAYS_IGNORE_ATTACKS)))

y=np.zeros(N_CHANS*max([len(_) for _ in y_sigs]))
for chan,y_ in enumerate(y_sigs):
    y[chan:len(y_)*N_CHANS:N_CHANS] = y_

y=common.normalize(y)*.99
y.tofile(OUT_FILE)
