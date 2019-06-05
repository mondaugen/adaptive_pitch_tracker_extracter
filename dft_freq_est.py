import numpy as np
from scipy import signal
import librosa

def pick_peaks(x,thresh=float('-inf')):
    dx=np.diff(x)
    dirx=np.zeros_like(dx)
    dirx[dx>=0]=1
    dirx[dx<0]=-1
    peaks=np.diff(dirx)
    peaks_i=np.where((peaks==-2)&(x[1:-1]>thresh))[0]+1
    return peaks_i

class peak_analyzer:
    def __init__(self,
        N_FFT,
        N_W,
        H=1,
        w_type='hann'):
        assert(N_FFT > N_W)
        start_i=(N_FFT-N_W)//2
        end_i=start_i+N_W
        self.w=signal.get_window(w_type,N_W)
        self.w/=np.sum(self.w)
        self.start_i=start_i
        self.end_i=end_i
        self.N_FFT=N_FFT
        self.N_W=N_FFT
        self.H=H
    def freqs_amps(self,x,max_n_peaks=20,peak_thresh=-80):
        """
        x is chunk of signal to analyze
        max_n_peaks is the maximum number of peaks to output (peaks are output
        starting at the highest amplitude)
        peak_thresh is a value in decibels, below which a peak is not considered
        returns initial phase, initial amplitude, instantaneous frequency,
        instantaneous amplitude decay
        """
        assert(len(x)==(self.N_W+self.H))
        xh0=np.zeros(self.N_FFT)
        xh1=np.zeros(self.N_FFT)
        xh0[self.start_i:self.end_i]=x[:self.N_W]*self.w
        xh1[self.start_i:self.end_i]=x[self.H:]*self.w
        Xh0=np.fft.rfft(xh0)
        Xh1=np.fft.rfft(xh1)
        Xh0_a=np.abs(Xh0_a)
        peaks_i=pick_peaks(Xh0_a,thresh=peak_thresh)
        peaks_i=peaks_i[np.argsort(Xh0_a[peaks_i])[::-1]][:max_n_peaks]
        Xh0_=Xh0[peaks_i]
        Xh1_=Xh1[peaks_i]
        ph=np.angle(Xh0_)
        a=np.abs(Xh0_)
        dph=(np.angle(Xh1_)-np.angle(Xh0_))/self.H
        da=np.power(np.abs(Xh1_/Xh0_),1./self.H)
        return (ph,a,dph,da)


