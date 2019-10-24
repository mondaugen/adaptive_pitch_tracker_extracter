import numpy as np
from scipy import signal
# Time stretch by estimated phase advancement using "ghost" windows.

def time_stretch_arb_times(
    x, # signal to time stretch
    t, # times at which to place windows
    H, # hop size
    W, # window size
    window_type='hann',
    synth_window_type=None):
    if synth_window_type is None:
        synth_window_type=window_type

    # check to see that all the times are valid
    # a valid time t is H <= t < len(x)
    # the first time can be 0 < t len(x) because it is just copied verbatim
    if not (np.all(t[1:] <= (len(x)-W))):
        raise IndexError('Some non-initial t greater than last index minus window length')
    if not ((t[0] >= 0) and (t[0] < (len(x)-W))):
        raise IndexError('Initial t not >= 0, or greater than last index minus window length')
    # pad beginning of signal with a hop size worth of 0s
    #x=np.concatenate((np.zeros(H,dtype=x.dtype),x))
    x=np.concatenate((np.random.standard_normal(H).astype(x.dtype)*1e-8,x))
    t_shift=H 
    # get analysis window
    w=signal.get_window(window_type,W)
    # get synthesis window
    sw=signal.get_window(synth_window_type,W)
    y_len=H*(len(t)-1)+W
    y=np.zeros(y_len)
    t0=t[0]+t_shift
    win_1=np.zeros_like(w)
    # calculate dividing out window, even for hops that aren't multiple of window size
    w_=np.concatenate((w,np.zeros(H-(W%H))))
    sw_=np.concatenate((sw,np.zeros(H-(W%H))))
    #win_div=np.sum(np.power(w,2).reshape((W//H,H)),axis=0)
    win_div=np.sum((w_*sw_).reshape((len(w_)//H,H)),axis=0)
    print('win_div:')
    print(win_div)
    # TODO: Why is the beginning so difficult?
    for h in np.arange(0,W,H):
        #win_1[:W-h]+=w[h:]
        _Y=np.fft.fft(np.concatenate((np.zeros(h),x[t0:t0+W-h]*w[h:])))
        if h == 0:
            Y_last=_Y.copy()
        y[:W-h]+=np.real(np.fft.ifft(Y_last))[h:]*sw[h:]
    y[:H]/=win_div
    h=H
    for t_ in t[1:]:
        t_=t_+t_shift
        X0=np.fft.fft(x[t_:t_+W]*w)
        X_H=np.fft.fft(x[t_-H:t_-H+W]*w)
        Y_last=X0*np.abs(X_H)/X_H*Y_last/np.abs(Y_last)
        y[h:h+W]+=np.real(np.fft.ifft(Y_last))*sw
        y[h:h+H]/=win_div
        h+=H

    return y
