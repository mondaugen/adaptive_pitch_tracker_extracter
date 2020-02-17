# Pitch-shift and or time-stretch multi-channel file by constant amounts

import psts
import common
import numpy as np
import attack_finder

IN_FILE=common.get_env('IN_FILE',default='/tmp/in.f64')
OUT_FILE=common.get_env('OUT_FILE',default='/tmp/out.f64')
SR=common.get_env('SR',default=16000,conv=int)
LMAX_FILT_RATE=common.get_env('LMAX_FILT_RATE',default=SR,conv=float)
N_CHANS=common.get_env('N_CHANS',default=1,conv=int)
W=common.get_env('W',default=1024,conv=int)
H=common.get_env('H',default=256,conv=int)
ATTACK_EST_MODE=common.get_env('ATTACK_EST_MODE',default='single-channel')
ATTACK_EST_CHAN=common.get_env('ATTACK_EST_CHAN',default=0,conv=int)
ATTACK_EST_W=common.get_env('ATTACK_EST_W',default=W,conv=int)
ATTACK_EST_H=common.get_env('ATTACK_EST_H',default=H,conv=int)

x=np.fromfile(IN_FILE,dtype='float64')

if ATTACK_EST_MODE == 'single-channel':
    # just obtain the attack times from the channel ATTACK_EST_CHAN and use them
    # for all channels
    x_attack_est=x[ATTACK_EST_CHAN::N_CHANS]
elif ATTACK_EST_MODE == 'all-channels':
    # obtain the attack times from the sum of all channels and use them
    # for all channels
    x_attack_est=np.zeros(len(x)//N_CHANS)
    for chan in range(N_CHANS):
        x_chan=x[chan::N_CHANS]
        x_attack_est[:len(x_chan)]+=x_chan
else:
    raise ValueError("Invalid ATTACK_EST_MODE %s" % (ATTACK_EST_MODE,))

attack_time_pairs=attack_finder.attacks_from_spectral_diff(
    x_attack_est,
    W=ATTACK_EST_W,
    H=ATTACK_EST_H,
    lmax_filt_rate=LMAX_FILT_RATE)
attack_times=np.array([b for a,b in attack_time_pairs])

y_sigs=[]

for chan in range(N_CHANS):
    y_sigs.append(psts.psts_const_amount(
        x[chan::N_CHANS],
        attack_times,
        SR=SR,
        W=W,
        H=H,
        ADSR_ATTACK=common.get_env('ADSR_ATTACK',default=0.1,conv=float),
        ADSR_RELEASE=common.get_env('ADSR_RELEASE',default=0.1,conv=float),
        TS=common.get_env('TS',default=1,conv=float),
        PS=common.get_env('PS',default=1,conv=float)))

y=np.zeros(N_CHANS*max([len(_) for _ in y_sigs]))
for chan,y_ in enumerate(y_sigs):
    y[chan:len(y_)*N_CHANS:N_CHANS] = y_

y.tofile(OUT_FILE)
