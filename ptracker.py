import numpy as np
from scipy.signal import lfilter, freqz
import matplotlib.pyplot as plt

class lpfilt:
    """
    Low pass filter.
    """
    def __init__(self,a):
        self.a = a
        self.b0 = (1-a)
        self.zi = None
    def _get_b_a(self):
        return ([self.b0],[1,self.a])#([self.b0,self.b0],[1,self.a])
    def filter(self,x):
        if self.zi is None:
            y,self.zi = lfilter(*self._get_b_a(),x,zi=[0])
        else:
            y,self.zi = lfilter(*self._get_b_a(),x,zi=self.zi)
        return y
    def plot_freqz(self,ax):
        w,h=freqz(*self._get_b_a(),whole=True)
        ax.plot(w,20*np.log10(np.abs(h)))

class dirfilt:
    """
    The "direction" filter.
    d = 1 -> positive
    d = -1 -> negative
    """
    def __init__(self,a,d=1):
        self.a = a
        self.arg_a = np.exp(complex("j")*np.pi*0.5*d)
        self.b0 = (1 - self.a)
        self.zi = None
    def _get_b_a(self):
        return ([self.b0],[1,self.a*self.arg_a])
    def filter(self,x):
        if self.zi is None:
            y,self.zi = lfilter(*self._get_b_a(),x,zi=[0])
        else:
            y,self.zi = lfilter(*self._get_b_a(),x,zi=self.zi)
        return y
    def plot_freqz(self,ax):
        w,h=freqz(*self._get_b_a(),whole=True)
        ax.plot(w,20*np.log10(np.abs(h)))

class apt:
    """
    Adaptive partial tracker.
    """
    def __init__(self,
        # Radius of low pass filter pole
        a_lp,
        # Radius of "positive" or "negative" filter pole
        a_pn,
        # Running power computation weighting coefficient
        alph_p,
        # instantaneous frequency correction weighting coefficient
        alph_w,
        # output amplitude smoothing
        alph_a=0.99
    ):
        self.lpf = lpfilt(a_lp)
        self.cpnf = dirfilt(a_pn,d=1)
        self.cnnf = dirfilt(a_pn,d=-1)
        self.alph_p = alph_p
        self.ppn = 0
        self.pnn = 0
        self.alph_w = alph_w
        self.alph_a = alph_a
        self.a_=None
    def proc(self,
        # Array of input values
        x,
        # initial guess 
        w0,
    ):
        """ Returns array a of amplitudes and array phi of phases for the complex exponential """
        g_n = complex('1')
        w = w0
        a = np.zeros((len(x),))
        phi = np.zeros((len(x),))
        corr = np.zeros((len(x),))
        err = np.zeros((len(x),),dtype='complex')
        for n,x_n in enumerate(x):
            e = x_n / g_n
            l = self.lpf.filter(np.array([e]))
            cpn = self.cpnf.filter(l)[0]
            cnn = self.cnnf.filter(l)[0]
            self.ppn = np.real(cpn*np.conj(cpn)*(1-self.alph_p)) + self.alph_p*self.ppn
            self.pnn = np.real(cnn*np.conj(cnn)*(1-self.alph_p)) + self.alph_p*self.pnn
            s = self.ppn - self.pnn
            w -= self.alph_w * s
            g_n *= np.exp(w*complex("j"))
            if self.a_ is None:
                self.a_ =  np.abs(l)
            else:
                self.a_ = np.abs(l)*(1-self.alph_a) + self.a_*self.alph_a
            a[n] = self.a_
            phi[n] = np.angle(g_n)
            corr[n] = self.alph_w * s
            err[n] = e
        return (a,phi,corr,err)
    def plot(self):
        fig,axs=plt.subplots(3,1)
        self.lpf.plot_freqz(axs[0])
        self.cpnf.plot_freqz(axs[1])
        self.cnnf.plot_freqz(axs[2])
        plt.show()

