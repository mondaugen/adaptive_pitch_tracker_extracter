# Time stretch by estimated phase advancement using "ghost" windows.

import numpy as np
from scipy import signal
import common 

def prepare_init_window(w_x_sw,H):
    """
    The first window is the product of the analysis window and the synthesis
    window, overlapped and added so that the initial frames (before the first
    frame) can just be added in without needing to do the whole overlap-and-add
    routine
    """
    W=len(w_x_sw)
    v=np.zeros_like(w_x_sw)
    for h in range(0,W,H):
        v[:W-h] += w_x_sw[h:]
    return v

def tsat_oob_fill_func_default(N):
    return np.random.standard_normal(N)*1e-8

def time_stretch_arb_times(
    x, # signal to time stretch
    t, # times at which to place windows
    H, # hop size
    W, # window size
    window_type='hann',
    synth_window_type=None,
    oob_fill_func=tsat_oob_fill_func_default):

    if synth_window_type is None:
        synth_window_type=window_type

    min_t=np.min(t)
    max_t=np.max(t)

    # pad signal so that all t's are valid
    # beginning of x is considered to be time 0
    t_off=-min(min_t,0)+H
    x=np.concatenate((
        oob_fill_func(t_off).astype(x.dtype),
        x,
        oob_fill_func(max(max_t + W - len(x),0)).astype(x.dtype)))

    # get analysis window
    w=signal.get_window(window_type,W)
    # get synthesis window
    sw=signal.get_window(synth_window_type,W)
    y_len=H*(len(t)-1)+W
    y=np.zeros(y_len)
    # Offset the times by the number of samples padded to the beginning of x
    t += t_off
    # calculate dividing out window, even for hops that aren't multiple of window size
    win_div=1/common.ola_shorten(w*sw,H)
    # calculate initial window so we don't need to do extra ola's to fill the initial buffer
    init_win=prepare_init_window(w*sw,H)
    # sum in "fake" inital overlap-and-add frames
    y[:W]=init_win*x[t[0]:t[0]+W]
    # store initial Y_last
    Y_last=np.fft.fft(x[t[0]:t[0]+W]*w)
    y[:H]*=win_div
    h=H
    for t_ in t[1:]:
        t_=t_
        X0=np.fft.fft(x[t_:t_+W]*w)
        X_H=np.fft.fft(x[t_-H:t_-H+W]*w)
        Y_last=X0*np.abs(X_H)/X_H*Y_last/np.abs(Y_last)
        y[h:h+W]+=np.real(np.fft.ifft(Y_last))*sw
        y[h:h+H]*=win_div
        h+=H

    return y
