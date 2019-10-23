import numpy as np
from scipy import signal
# Time stretch by estimated phase advancement using "ghost" windows.

def time_stretch_arb_times(
    x, # signal to time stretch
    t, # times at which to place windows
    H, # hop size
    W, # window size
    window_type='hann'):

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
    w=signal.get_window(window_type,W)
    y_len=H*(len(t)-1)+W
    y=np.zeros(y_len)
    t0=t[0]+t_shift
    win_1=np.zeros_like(w)
    s=np.zeros_like(y)
    # TODO: Why is the beginning so difficult?
    for h in np.arange(0,W,H):
        #win_1[:W-h]+=w[h:]
        _Y=np.fft.fft(np.concatenate((np.zeros(h),x[t0:t0+W-h]*w[h:])))
        if h == 0:
            Y_last=_Y.copy()
        y[:W-h]+=np.real(np.fft.ifft(Y_last))[h:]*w[h:]
        s[:W-h]+=np.power(w[h:],2)

    h=H
    for t_ in t[1:]:
        t_=t_+t_shift
        X0=np.fft.fft(x[t_:t_+W]*w)
        X_H=np.fft.fft(x[t_-H:t_-H+W]*w)
        Y_last=X0*np.abs(X_H)/X_H*Y_last/np.abs(Y_last)
        y[h:h+W]+=np.real(np.fft.ifft(Y_last))*w
        s[h:h+W]+=np.power(w,2)
        h+=H
    y/=s

    return y
        
