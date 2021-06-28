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

def fractional_get_window_dft(wintype,N_win,N_dft,oversample=1):
    """
    Get the Fourier transform of length N_dft of a window of length N_win of
    type wintype.  N_win can be fractional: the Fourier transforms of the
    greatest smaller and smallest greater windows of even length are found and
    interpolated linearly.
    oversample is the number of times the fourier transform is oversampled
    """
    N_win_f=2*floor(N_win/2)
    N_win_c=2*ceil(N_win/2)
    if N_win_f == N_win_c:
        frac = 0
    else:
        frac=(N_win - N_win_f)/(N_win_c - N_win_f)
    w0=get_padded_window(wintype,N_win_f,N_dft)
    w1=get_padded_window(wintype,N_win_c,N_dft)
    bins=np.arange(N_dft*oversample)/oversample
    F=np.exp(-j*2*np.pi*bins[:,None]*np.arange(N_dft)/N_dft)/N_dft
    W0=F@w0
    W1=F@w1
    return W0 + frac * (W1 - W0)

# This should be compatible with dftX, dft, dft_dv accepting R argument
class multi_q_window:
    def __init__(
    # Maximum window length (length of DFT)
    N_max,
    # Desired Q for specific, unique v (v normalized frequency), ascending in v and
    # as tuples e.g., [(v0,Q0),(v1,Q1),....]
    Q_v,
    # The lobe radius in bins (corresponding to window type)
    B_l,
    # The window type (e.g., 'blackman')
    wintype,
    # The number of times the fourier transform of the windows are oversampled
    oversample=16):
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
        self.N_max=N_max
        self.Q_v=Q_v
        self.B_l=B_l
        self.wintype=wintype
        self.oversample=oversample
        # The minimum Q possible at a given v for N_max is
        Q_min=lambda v: (v + B_l / N_max) / v
        # Check for overzealous Q
        if any([Q_min(v) > Q for v,Q in Q_v]):
            raise ValueError('Infeasible Q(v) for given N_max.')
        self.v=np.array([v for v,Q in Q_v])
        self.Q=np.array([Q for v,Q in Q_v])
        # Ideal window length (can be fractional)
        self.N_v=B_l/(self.v*(self.Q - 1))
        # Find the DFTs of windows for each N_v
        self.Ws=np.vstack([fractional_get_window_dft(wintype,n_v,N_max,oversample=oversample)
                      for n_v in self.N_v])
        # Find the new lobe radii caused by the smaller N_v
        self.B_l_v=B_l*N_max/self.N_v
        self.win_dft_interp=interpolate.interp2d(np.arange(N_max),self.v,self.Ws)
        self.lobe_radius_interp=(lambda v_: np.ceil(
            interpolate.interp1d(self.v,self.B_l_v)(v_)))
    def __len__(self):
        return self.N_max
    def R(self,v):
        """
        Get a sparse matrix R that when right-multiplied by the fourier
        transform of a signal of length N_win gives the values of the fourier
        transform of the signal at the bins v.
        v are the normalized frequencies.
        """
        R=sparse.lil_matrix((len(v),self.N_win),dtype='complex128')
        b=v*self.N_win
        b_rounded=np.round(b).astype('int')
        b_frac=b_rounded-b
        # Figure out what bins are worth looking up
        lobe_radii=self.lobe_radius_interp(v)
        lu_bins=[np.arange(-lr,lr+1) for lr in lobe_radii]
        b_ranges=[lub - bf for bf,lub in zip(b_frac,lu_bins)]
        b_indices=[(lub - br) % self.N_max for br,lub in zip(b_rounded,lu_bins)]
        r_indices=[idx * np.ones(len(lub)) for idx,lub in enumerate(lu_bins)]
        vs=[v_ * np.ones(len(lub)) for v_,lub in zip(v,lu_bins)]
        R[np.concatenate(r_indices),np.concatenate(b_indices)]=self.win_dft_interp(
            np.concatenate(b_ranges),np.concatenate(vs))
        return R.tocsr()

