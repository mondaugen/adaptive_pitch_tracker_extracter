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
        self.N_W=N_W
        self.H=H

    def dft(self,x):
        assert(len(x)==(self.N_W))
        xh0=np.zeros(self.N_FFT)
        xh0[self.start_i:self.end_i]=x[:self.N_W]*self.w
        Xh0=np.fft.rfft(xh0)
        return Xh0

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
        Xh0_a=np.abs(Xh0)
        peaks_i=pick_peaks(Xh0_a,thresh=np.power(10,peak_thresh/20))
        peaks_i=peaks_i[np.argsort(Xh0_a[peaks_i])[::-1]][:max_n_peaks]
        Xh0_=Xh0[peaks_i]
        Xh1_=Xh1[peaks_i]
        ph=np.angle(Xh0_)
        ph1=np.angle(Xh1_)
        a=np.abs(Xh0_)
        # phase unwrap calc
        # use bin frequencies of peaks to get an estimate of the frequency
        w_k=peaks_i/self.N_FFT*2*np.pi
        # we now see how many times we need to unwrap in order to minimize the
        # difference between the estimated frequency and the bin frequency
        A=(self.H*w_k - ph1 + ph)/self.H
        # how many times to unwrap by 2*pi
        M=np.round(self.H*A/(2*np.pi))
        dph=(ph1+M*np.pi*2-ph)/self.H
        da=np.power(np.abs(Xh1_/Xh0_),1./self.H)
        return (ph,a,dph,da)

    def nearest_peaks(self,x,w_ranges,default_ph=0,default_a=1e-5):
        """
        Output only peaks that fall in one of the ranges in w_ranges.
        If no peak falls in a range, put the mid-point of the range as the
        frequency, default_ph and default_a for the phase and amplitude, and 1
        for the decay rate.
        Returns ph_,a_,dph_,da_ for the peaks
        """
        ph,a,dph,da=self.freqs_amps(x,max_n_peaks=self.N_FFT,peak_thresh=20*np.log10(default_a))
        w_midpoints=np.mean(w_ranges,axis=1)
        diffs=np.abs(np.subtract.outer(w_midpoints,dph))
        ord_diffs=np.argsort(diffs)
        ph_=np.ones(w_midpoints.shape)*default_ph
        a_=np.ones(w_midpoints.shape)*default_a
        dph_=w_midpoints
        da_=np.ones(w_midpoints.shape)
        found_peaks=(dph[ord_diffs[:,0]]<=w_ranges[:,1])&(dph[ord_diffs[:,0]]>=w_ranges[:,0])
        orig_data_i=ord_diffs[:,0][found_peaks]
        ph_[found_peaks]=ph[orig_data_i]
        a_[found_peaks]=a[orig_data_i]
        dph_[found_peaks]=dph[orig_data_i]
        da_[found_peaks]=da[orig_data_i]
        return (ph_,a_,dph_,da_)

    def synth_peaks(self,ph,a,dph,da,ignore_da=True,return_th_a=False):
        """
        Using the pairs of (ph,a,dph,da) in each field, synthesize a segment of a sinusoid.
        if ignore_da is true, amplitude is linearly interpolated between a,
        otherwise da is taken into account when designing an interpolating
        polynomial for the amplitude.
        """
        if ignore_da==False:
            raise NotImplementedError
        else:
            _da=np.power(a[:,1]/(a[:,0]+1e-8),1./self.H)
        # Find multiple of 2pi that unwraps the phase the most smoothly
        M=np.round((ph[:,0]+dph[:,0]*self.H-ph[:,1]+self.H*0.5*(dph[:,1]-dph[:,0]))/(np.pi*2))
        # compute interpolating polynomial coefficients
        alph=(3/(self.H*self.H)*(ph[:,1]-ph[:,0]-dph[:,0]*self.H+2*np.pi*M)
            - 1/self.H*(dph[:,1]-dph[:,0]))
        beta=((dph[:,1]-dph[:,0])/(self.H*self.H) 
            - 2/(self.H*self.H*self.H)*(ph[:,1]-ph[:,0]-dph[:,0]*self.H+2*np.pi*M))
        # evaluate phase function
        n=np.multiply.outer(np.ones(len(alph)),np.arange(self.H))
        th=beta[:,None]*n
        th=n*(alph[:,None]+th)
        th=n*(dph[:,0][:,None]+th)
        th+=ph[:,0][:,None]
        # evaluate amplitude function
        amp=a[:,0][:,None]*np.power(_da[:,None],n)
        if return_th_a:
            return (amp,th)
        return np.sum(amp*np.cos(th),axis=0)

    def synth_peaks_trans_stretch(self,ph,a,dph,da,ignore_da=True,ptrans=1,tstretch=1):
        """
        Using the pairs of (ph,a,dph,da) in each field, synthesize a segment of a sinusoid.
        if ignore_da is true, amplitude is linearly interpolated between a,
        otherwise da is taken into account when designing an interpolating
        polynomial for the amplitude.
        Here you can change the pitch by rescaling n when synthesizing the
        phase, or change the time by rescaling the hop size.

        The way pitch transposition and or time-stretching works is by rescaling
        or extending n.
        For pitch transposition, when synthesizing the phase, we rescale n so
        that the phase advances more or less than it would with no
        transposition. Then we also return the final phase so that the next
        synthesis step can use this corrected phase. The amplitude envelope remains the same.
        For time stretching, we extend or contract n, depending on whether we
        are stretching or contracting time. we then synthesize into this
        differently sized buffer. The altered phase is returned, too. The
        amplitude envelope is contracted or expanded accordingly as well.
        We also return the error in contracting / expanding the buffer size
        because this will be integer valued. When enough errors accumulate, this
        can be used to apply the correct amount of time-stretch in a later
        frame.
        """
        if ignore_da==False:
            raise NotImplementedError
        else:
            _da=np.power(a[:,1]/(a[:,0]+1e-8),1./self.H)
        # Find multiple of 2pi that unwraps the phase the most smoothly
        M=np.round((ph[:,0]+dph[:,0]*self.H-ph[:,1]+self.H*0.5*(dph[:,1]-dph[:,0]))/(np.pi*2))
        # compute interpolating polynomial coefficients
        alph=(3/(self.H*self.H)*(ph[:,1]-ph[:,0]-dph[:,0]*self.H+2*np.pi*M)
            - 1/self.H*(dph[:,1]-dph[:,0]))
        beta=((dph[:,1]-dph[:,0])/(self.H*self.H) 
            - 2/(self.H*self.H*self.H)*(ph[:,1]-ph[:,0]-dph[:,0]*self.H+2*np.pi*M))
        # evaluate phase function
        N_=np.round(self.H*tstretch)
        N_err=self.H*tstretch-N_
        # NOTE: This tstretch_err is not useful
        tstretch_err=self.H*tstretch/N_
        n_a=np.arange(N_+1)/N_*self.H
        n_ph=np.arange(N_+1)*ptrans
        n=np.multiply.outer(np.ones(len(alph)),n_ph)
        na=np.multiply.outer(np.ones(len(alph)),n_a)
        th=beta[:,None]*n
        th=n*(alph[:,None]+th)
        th=n*(dph[:,0][:,None]+th)
        th+=ph[:,0][:,None]
        # evaluate amplitude function
        amp=a[:,0][:,None]*np.power(_da[:,None],na)
        th_ret=th[:,-1]
        #while np.any(th_ret>=(2*np.pi)) or np.any(th_ret<(2*np.pi)):
        #    th_ret[th_ret>=2*np.pi]-=2*np.pi
        #    th_ret[th_ret<2*np.pi]+=2*np.pi
        th_ret = th_ret % (2*np.pi)
        return (np.sum(amp*np.cos(th),axis=0)[:-1],amp[:,-1],th_ret,tstretch_err,amp,th)

    def synth_peaks_trans(self,ph,a,dph,da,ptrans=1,ignore_da=True):
        if ignore_da==False:
            raise NotImplementedError
        else:
            _da=np.power(a[:,1]/(a[:,0]+1e-8),1./self.H)
        # Adjust the predicted phase using the estimated frequency
        ph_adj=0.5*(1-ptrans)*(dph[:,0]+dph[:,1])*self.H
        #print(ph_adj)
        ph[:,1]-=ph_adj
        dph*=ptrans
        # Find multiple of 2pi that unwraps the phase the most smoothly
        M=np.round((ph[:,0]+dph[:,0]*self.H-ph[:,1]+self.H*0.5*(dph[:,1]-dph[:,0]))/(np.pi*2))
        # compute interpolating polynomial coefficients
        alph=(3/(self.H*self.H)*(ph[:,1]-ph[:,0]-dph[:,0]*self.H+2*np.pi*M)
            - 1/self.H*(dph[:,1]-dph[:,0]))
        beta=((dph[:,1]-dph[:,0])/(self.H*self.H) 
            - 2/(self.H*self.H*self.H)*(ph[:,1]-ph[:,0]-dph[:,0]*self.H+2*np.pi*M))
        # evaluate phase function
        n=np.multiply.outer(np.ones(len(alph)),np.arange(self.H+1))
        th=beta[:,None]*n
        th=n*(alph[:,None]+th)
        th=n*(dph[:,0][:,None]+th)
        th+=ph[:,0][:,None]
        # evaluate amplitude function
        amp=a[:,0][:,None]*np.power(_da[:,None],n)
        return (np.sum(amp*np.cos(th),axis=0)[:-1],th[:,-1],amp[:,:-1],th[:,:-1])
