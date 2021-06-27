# calculate windows have Q specified for specific frequencies

import numpy as np
from scipy import interpolate, signal
from common import is_pow_2
from math import floor, ceil

def get_padded_window(wintype,N,N_max):
    # Only works if N_max is power of 2 and N is even
    assert is_pow_2(N_max)
    assert N <= N_max
    assert (N % 2) == 0
    w=np.zeros(N_max)
    w_=signal.get_window(wintype,N)
    w[:N//2]=w_[N//2:]
    w[-N//2:]=w_[:N//2]
    return w

def fractional_get_window_dft(wintype,N_win,N_dft):
    """
    Get the Fourier transform of length N_dft of a window of length N_win of
    type wintype.  N_win can be fractional: the Fourier transforms of the
    greatest smaller and smallest greater windows of even length are found and
    interpolated linearly.
    """
    N_win_f=2*floor(N_win/2)
    N_win_c=2*ceil(N_win/2)
    if N_win_f == N_win_c:
        frac = 0
    else:
        frac=(N_win - N_win_f)/(N_win_c - N_win_f)
    w0=get_padded_window(wintype,N_win_f,N_dft)
    w1=get_padded_window(wintype,N_win_c,N_dft)
    W0=np.fft.fft(w0)
    W1=np.fft.fft(w1)
    return W0 + frac * (W1 - W0)

def multi_q_window(
# Maximum window length (length of DFT)
N_max,
# Desired Q for specific, unique v (v normalized frequency), ascending in v and
# as tuples e.g., [(v0,Q0),(v1,Q1),....]
Q_v,
# The lobe radius in bins (corresponding to window type)
B_l,
# The window type (e.g., 'blackman')
wintype):
    """
    Returns 2 interpolators.
    The first is an instance of interpolate.interp2d that can be evaluated with x
    and y, where x is the bin number in [0,N_max) and y is the normalized
    frequency in (-inf,+inf). Out of bounds y are made to be the minimum or
    maximum frequency speicifed.
    The second is a function that can be evalutated with the normalized
    frequency (accepts np.arrays). This gives the approximate lobe radius in
    bins, rounded up.
    """
    # The minimum Q possible at a given v for N_max is
    Q_min=lambda v: (v + B_l / N_max) / v
    # Check for overzealous Q
    if any([Q_min(v) > Q for v,Q in Q_v]):
        raise ValueError('Infeasible Q(v) for given N_max.')
    v=np.array([v for v,Q in Q_v])
    Q=np.array([Q for v,Q in Q_v])
    # Ideal window length (can be fractional)
    N_v=B_l/(v*(Q - 1))
    # Find the DFTs of windows for each N_v
    Ws=np.vstack([fractional_get_window_dft(wintype,n_v,N_max) for n_v in N_v])
    # Find the new lobe radii caused by the smaller N_v
    B_l_v=B_l*N_max/N_v
    win_dft_interp=interpolate.interp2d(np.arange(N_max),v,Ws)
    lobe_radius_interp=lambda v_: np.ceil(interpolate.interp1d(v,B_l_v)(v_))
    return (win_dft_interp,lobe_radius_interp)
