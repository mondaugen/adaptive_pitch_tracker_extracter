# Apply windows in the frequency domain
# This is possible because of Parseval's theorem.
# This is useful when it is desired to compute the DFT of bins that are not at
# integer multiples of the window's frequency.

import numpy as np
from scipy import signal, interpolate, sparse

import pdb

j=complex('j')

def calc_window_with_radius(window_name,lobe_radius):
    def f(N_win,oversamp):
        l=N_win
        w=signal.get_window(window_name,N_win)
        w/=np.sum(w)
        evalradius=lobe_radius*oversamp
        evalbins=np.arange(-evalradius,evalradius+1)
        bins=evalbins/oversamp
        # Oversampled fourier transform matrix
        # Divide by N_win to remove the fourier transform constant
        F=np.exp(-j*2*np.pi*bins[:,None]*np.arange(N_win)/N_win)/N_win
        vals=(F@w)[:len(evalbins)]
        return (bins,vals)
    return f

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
        win_type is the type of window, e.g., 'blackman'
        interpolator is the interpolation method (currently only 'linear' is
        officially supported)
        oversample is the number of times the DFT of the window is oversampled
        when computing the values for the interpolator.
        """
        self.N_win=N_win
        bins,vals=window_types[win_type]['calc'](N_win,oversample)
        self.bins=bins
        self.vals=vals
        self.W_lookup=interpolate.interp1d(self.bins,self.vals,kind=interpolator,
        fill_value=0.,bounds_error=False)
        lobe_radius=window_types[win_type]['lobe_radius']
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
        b_ranges=-b_frac[:,None]+self.lu_bins
        b_indices=(-b_rounded[:,None]+self.lu_bins).astype('int') % self.N_win
        r_indices=np.arange(len(v))[:,None]*np.ones(len(self.lu_bins))
        R[r_indices.flatten(),b_indices.flatten()]=self.W_lookup(b_ranges.flatten())
        return R.tocsr()

def calc_X(x):
    X=np.fft.fft(x)
    return X

def calc_dX(x):
    N=len(x)
    x_=x*j*2*np.pi*np.arange(N)
    dX=calc_X(x_)
    return dX
        
def dftX(X,R):
    RX=R@np.conj(X)
    return RX

def dft(x,R):
    X=calc_X(x)
    return dftX(X,R)

def dft_dv(x,R):
    x_=calc_dX(x)
    return dft(x_,R)

    
