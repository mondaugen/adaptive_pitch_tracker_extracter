# Apply windows in the frequency domain
# This is possible because of Parseval's theorem.
# This is useful when it is desired to compute the DFT of bins that are not at
# integer multiples of the window's frequency.

import numpy as np
from scipy import signal, interpolate, sparse
from common import is_pow_2
from math import floor, ceil
from functools import partial
# TODO: some_sig, dftdk shouldn't be imported from a test folder, currently
# requires test/gradient_partial_tracker/ in PYTHONPATH
from some_sig import multiply_ramp
from dftdk import half_shift_x

j=complex('j')

def calc_window_with_radius(window_name,lobe_radius):
    def f(N_win,oversamp):
        W=fractional_get_window_dft(window_name,N_win,N_win,oversample=oversamp,real=False,symmetrical_padding=False)
        # TODO: real=False here as we can't interpolate complex numbers as
        # easily as real numbers, and there should be no real reason to use
        # complex numbers (we can only support symmetric windows, which is good
        # enough)
        evalradius=lobe_radius*oversamp
        evalbins=np.arange(-evalradius,evalradius+1)
        bins=evalbins/oversamp
        # Oversampled fourier transform matrix
        # Divide by N_win to remove the fourier transform constant
        vals=np.concatenate((W[-evalradius:],W[:evalradius+1]))
        return (bins,vals)
    return f

calc_blackmanharris=calc_window_with_radius('blackmanharris',4)
calc_blackman=calc_window_with_radius('blackman',3)
calc_hann=calc_window_with_radius('hann',2)

window_types = {
    # The name of the window type
    'blackman' : {
        # The distance to the first zero after the center of the mainlobe in bins.
        'lobe_radius' : 3,
        # A function taking the arguments 'N_win' and 'oversamp'
        # 'N_win' is the length of the window in samples
        # 'oversamp' is the number of times of oversampling, so oversamp=2 would
        # give values halfway between the bins as well as the bin values
        # this function must return bins and vals, where bins are the bins at
        # which the DFT of the window was evaluated and vals the values there
        # NOTE: This function should return real values
        'calc' : calc_blackman
    },
    'blackmanharris' : {
        'lobe_radius' : 4,
        'calc' : calc_blackmanharris
    },
    'hann' : {
        'lobe_radius': 2,
        'calc': calc_hann
    }
}

class freq_dom_window:
    def __init__(self,N_win,win_type,oversample,interpolator='linear'):
        """
        Make a freq_dom_window instance.
        N_win is the length of the window in samples.
        win_type is a string giving the type of window, e.g., 'blackman', or can
        be a callable object
        accepting N_win (the length of the window) and oversamp (the
        oversampling factor: N_win*oversamp gives the length the window is
        padded to before performing the transform). Calling this object returns
        at tuple that are the bins that were evaluated in the DFT and the values
        of those bins. The peak value of the values should be 1 (typically) so
        as not to change the gain of observered sinusoids. It also must give the
        lobe_radius when looking up win_type.lobe_radius.
        interpolator is the interpolation method (currently only 'linear' is
        officially supported)
        oversample is the number of times the DFT of the window is oversampled
        when computing the values for the interpolator.
        """
        self.N_win=N_win
        if '__call__' in dir(win_type):
            bins,vals=win_type(N_win,oversample)
            lobe_radius=win_type.lobe_radius
        else:
            bins,vals=window_types[win_type]['calc'](N_win,oversample)
            lobe_radius=window_types[win_type]['lobe_radius']
        self.bins=bins
        # Divide by N_win to cancel out fourier transform constant
        self.vals=vals/self.N_win
        self.W_lookup=interpolate.interp1d(self.bins,self.vals,kind=interpolator,
        fill_value=0.,bounds_error=False)
        self.lu_bins=np.arange(-lobe_radius,lobe_radius+1)
    def __len__(self):
        return self.N_win
    def R(self,v):
        """
        Get a sparse matrix R that when right-multiplied by the fourier
        transform of a signal of length N_win gives the values of the fourier
       transform of the signal at the bins v.
        v are the normalized frequencies.
        """
        R=sparse.lil_matrix((len(v),self.N_win),dtype='complex128')
        b=v*self.N_win
        b_rounded=np.round(b)
        b_frac=b_rounded-b
        b_ranges=b_frac[:,None]+self.lu_bins
        b_indices=(b_rounded[:,None]+self.lu_bins).astype('int') % self.N_win
        r_indices=np.arange(len(v))[:,None]*np.ones(len(self.lu_bins))
        R[r_indices.flatten(),b_indices.flatten()]=self.W_lookup(b_ranges.flatten())
        return R.tocsr()

def calc_X(x,half_shift=True):
    N=len(x)
    if half_shift:
        _x=half_shift_x(x,N)
    else:
        _x=x
    X=np.fft.fft(_x)
    return X

def calc_dX(x):
    N=len(x)
    _x=half_shift_x(x,N)
    dX=multiply_ramp(_x,N,1)
    return dX
        
def dftX(X,R):
    """ R.shape[1] must equal len(X) """
    RX=np.conj(R)@X
    return RX

def dft(x,R,half_shift=True):
    """ R.shape[1] must equal len(x) """
    X=calc_X(x,half_shift=half_shift)
    return dftX(X,R)

def dft_dv(x,R):
    """ R.shape[1] must equal len(x) """
    dX=calc_dX(x)
    return dft(dX,R,half_shift=False)

    

def get_padded_window(wintype,N,N_max,symmetrical=True):
    # Only works if N_max is power of 2 and N is even
    assert is_pow_2(N_max)
    assert N <= N_max
    assert (N % 2) == 0
    w=np.zeros(N_max)
    w_=signal.get_window(wintype,N)
    if symmetrical:
        w[:N//2]=w_[N//2:]
        w[-N//2:]=w_[:N//2]
    else:
        w[:N]=w_
    return w

def fractional_get_window_dft(wintype,N_win,N_dft,oversample=1,real=True,symmetrical_padding=True):
    """
    Get the Fourier transform of length N_dft of a window of length N_win of
    type wintype.  N_win can be fractional: the Fourier transforms of the
    greatest smaller and smallest greater windows of even length are found and
    interpolated linearly.
    oversample is the number of times the fourier transform is oversampled
    If real is True, output the real part only, otherwise output the complex spectrum.
    If symmetrical_padding is True and wintype is symmetrical (usual case), then
    the windows are zero padded such that the "periodized" window function is an
    even function, and the spectrum of the window is purely real and setting
    real=True discards no information.
    """
    N_dft *= oversample
    N_win_f=2*floor(N_win/2)
    N_win_c=2*ceil(N_win/2)
    if N_win_f == N_win_c:
        frac = 0
    else:
        frac=(N_win - N_win_f)/(N_win_c - N_win_f)
    w0=get_padded_window(wintype,N_win_f,N_dft,symmetrical=symmetrical_padding)
    w1=get_padded_window(wintype,N_win_c,N_dft,symmetrical=symmetrical_padding)
    W0=np.fft.fft(w0)/np.sum(w0)
    W1=np.fft.fft(w1)/np.sum(w1)
    ret = W0 + frac * (W1 - W0)
    if real:
        ret = np.real(ret)
    return ret

def calc_Q_min(v,B_l,N):
    return (v + B_l / N) / v

# calculate windows that have Q specified for specific frequencies
# This should be compatible with dftX, dft, dft_dv accepting R argument
class multi_q_window:
    def __init__(
    self,
    # Maximum window length (length of DFT)
    N_max,
    # Lobe radius multiplier for specific unique v (v normalized frequency),
    # ascending in v and as tuples e.g., [(v0,m0),(v1,m1),....]
    # TODO: Currently need at least 4 points for interpolation, could we use
    # lower order splines?
    m_v,
    # The window type (e.g., 'blackman')
    wintype,
    # The number of times the fourier transform of the windows are oversampled
    oversample=16,
    # The interpolation order on the grid of windows
    v_q_interp_order=3):
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
        B_l=window_types[wintype]['lobe_radius']
        self.N_max=N_max
        self.m_v=m_v
        self.B_l=B_l
        self.wintype=wintype
        self.oversample=oversample
        self.N_os=self.N_max*self.oversample
        assert all(0 <= v <= 1 for v,_ in m_v)
        assert all(m > 0 for _,m in m_v)
        self.v=np.array([v for v,m in m_v])
        self.m=np.array([m for v,m in m_v])
        # Ideal window length (can be fractional)
        self.N_v=N_max/self.m
        # Find the DFTs of windows for each N_v
        self.Ws=np.vstack([fractional_get_window_dft(wintype,n_v,N_max,oversample=oversample)
                      for n_v in self.N_v])
        # Find the new lobe radii caused by the smaller N_v
        self.B_l_v=B_l*self.m
        # TODO: What you want is interpolate.RectBivariateSpline ... or is it?
        self.win_dft_interp=interpolate.RectBivariateSpline(
        self.v,
        np.arange(self.N_os)/oversample - N_max//2,
        np.hstack((self.Ws[:,-self.Ws.shape[1]//2:],
                   self.Ws[:,:self.Ws.shape[1]//2])),
        kx=v_q_interp_order,
        ky=v_q_interp_order)
        self.lobe_radius_interp=(lambda v_: np.ceil(
            interpolate.interp1d(self.v,self.B_l_v)(v_)))
    def __len__(self):
        return self.N_max
    def R(self,v):
        """
        Get a sparse matrix R that when right-multiplied by the fourier
        transform of a signal of length N_max gives the values of the fourier
        transform of the signal at the bins v.
        v are the normalized frequencies.
        """
        R=sparse.lil_matrix((len(v),self.N_max),dtype='complex128')
        b=v*self.N_max
        b_rounded=np.round(b).astype('int')
        b_frac=b_rounded-b
        # Figure out what bins are worth looking up
        lobe_radii=self.lobe_radius_interp(v)
        lu_bins=[np.arange(-lr,lr+1) for lr in lobe_radii]
        b_ranges=[lub + bf for bf,lub in zip(b_frac,lu_bins)]
        b_indices=[(self.N_max - (lub + br)) % self.N_max for br,lub in zip(b_rounded,lu_bins)]
        r_indices=[idx * np.ones(len(lub)) for idx,lub in enumerate(lu_bins)]
        vs=[v_ * np.ones(len(lub)) for v_,lub in zip(v,lu_bins)]
        a=-self.N_max/2
        fdt_shift=[np.exp(-2*np.pi*j*a*(b_+lub)/self.N_max) for b_,lub in zip(b,lu_bins)]
        # TODO: the values getting written into R need to be multiplied by a
        # complex exponential sequence as to shift them in time (the window they
        # represent the fourier transform of is not one at the correct time)
        # specificaly I think the windows need to be shifted back half a DFT
        # length (the window length multiplied by oversample) in time.
        R[np.concatenate(r_indices),np.concatenate(b_indices)]=np.conj(self.win_dft_interp(
            np.concatenate(vs),np.concatenate(b_ranges),grid=False)*np.concatenate(fdt_shift))
        return R.tocsr()

# TODO: First try getting a sort of constant Q transform of an STFT: the
# difference is that the bins are still centred on the normalized frequencies
# 0/N, 1/N, ... , (N-1)/N rather than arbitrary frequencies (often an equal
# tempered scale or similar). This should widen the main lobes of the higher
# frequencies.
