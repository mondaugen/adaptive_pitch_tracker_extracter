import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def default_nharms(freq):
    return int(np.floor(0.5/freq))

def default_harm(freq,n):
    return freq*(n+1)

def default_weight(freq,n):
    return np.ones_like(n)

def default_variance(freq,n):
    # variance is chosen so that the gaussian has a value around 1% of its peak half way between 2 harmonics
    return np.ones_like(n)*(freq/(2*3))**2

def default_window(N):
    return signal.get_window('hann',N,fftbins=True)

def pick_peaks(x,thresh=float('-inf')):
    dx=np.diff(x,axis=0)
    dirx=np.zeros_like(dx)
    dirx[dx>=0]=1
    dirx[dx<0]=-1
    peaks=np.diff(dirx,axis=0)
    xi,yi=np.where((peaks==-2)&(x[1:-1,:]>thresh))
    iord=np.argsort(yi)
    peaks_i=(xi[iord]+1,yi[iord])
    return peaks_i

def peak_matrix(x,thresh=float('-inf')):
    """ Takes x and returns a matrix of the same size where its entries are 1 if
    the row index corresponds to a local maximum in a particular column,
    """
    dx=np.diff(x,axis=0)
    dirx=np.zeros_like(dx)
    dirx[dx>=0]=1
    dirx[dx<0]=-1
    peaks=np.diff(dirx,axis=0)
    res=np.zeros_like(x)
    res[1:-1,:][(peaks==-2)&(x[1:-1,:]>thresh)]=1
    return res

class rab_pitch:
    """
    Estimate pitches using Rabiner's technique.
    Pick peaks in the spectrum, create a sum of delta functions offset by the
    peak frequencies, then find the inner product of this and a sum of gaussian
    functions offset by harmonics of some fundamental.
    """
    def __init__(self,
        # the set of frequencies
        # in 2pi radians/s (i.e., 0.5 is nyquist)
        freqs,
        # function given the frequency that gives the number of harmonics
        # must return integer
        # e.g., lambda freq: int(np.floor(0.5/freq))
        nharms=default_nharms,
        # a function passed the harmonic number and fundamental frequency and
        # giving the frequency of the harmonic
        # e.g., lambda freq, n : freq*(n+1)
        # this function is passed in values from 0 to nharms-1
        harm=default_harm,
        # a function which is passed the harmonic number and gives the
        # weight of the gaussian
        # e.g., lambda n : 1/(n+1)
        # n will range from 0 to the frequency such that freq*n < 0.5
        weight=default_weight, 
        # a function passed the harmonic number which gives the variance (width)
        # of the gaussian
        # e.g. lambda n : 1/(n+1)
        # but this seems more arbitrary
        # n will range from 0 to the frequency such that freq*n < 0.5
        variance=default_variance,
        # the number of points in the created look-up table
        ntabpoints=4096,
        # the frame lengths
        N_frame=4096,
        # the window length (must be <= N_frame)
        # it is zero padded so it works.
        N_window=2048,
        # see default_window for example
        window=default_window,
        # the threshold in dB over which a peak is considered
        peak_thresh=float('-inf')):

        self.tables=np.ndarray((len(freqs),ntabpoints))
        pts=np.arange(0,1,1/ntabpoints)*0.5
        for tabn,f in enumerate(freqs):
            nharm=np.arange(nharms(f))
            harms=harm(f,nharm)
            weights=weight(f,nharm)
            variances=variance(f,nharm)
            x=np.exp(
                -np.power(
                    np.add.outer(-harms,pts),
                    2)/(2*variances[:,None]))*weights[:,None]
            x=np.sum(x,axis=0)
            self.tables[tabn,:]=x
        # interpolator using the tables to look up values
        self.inter=interp1d(
        pts,
        self.tables,
        bounds_error=False,
        fill_value=0)
        # window for dft
        self.w=np.zeros((N_frame))
        w_start=(N_frame-N_window)//2
        w_end=N_frame-w_start
        self.w[w_start:w_end]=window(N_window)
        self.w/=np.sum(self.w)
        self.peak_thresh=10**(peak_thresh/20)
        self.N_frame=N_frame
    def _plot_tables(self):
        plt.plot(np.multiply.outer(np.ones(self.tables.shape[0]),
            np.arange(self.tables.shape[1])).T,self.tables.T)
    def __call__(self,x):
        """
        x must be arranged into frames of length len(self.w)
        The frames must be the columns of x.
        Returns a matrix with each row indicating a particular note score and
        each column an fft frame.
        """
        X=np.abs(np.fft.rfft(x*self.w[:,None],axis=0))
        peakmat=peak_matrix(X,thresh=self.peak_thresh)
        print(peakmat.shape)
        # Create a matrix where its values are the normalized frequency of the
        # peak if a peak exists or otherwise -1
        mones=peakmat*np.arange(X.shape[0])[:,None]/self.N_frame
        mones[mones<=0]=-1
        # transpose and use to look up interpolated peak values
        note_scores=np.sum(self.inter(mones.T)*((X*peakmat).T[None,:]),axis=2)
        return note_scores
