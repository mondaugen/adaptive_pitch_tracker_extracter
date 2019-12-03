# LFO stuff
import numpy as np
from scipy import signal

class chirp:

    def __init__(self,
        # length of signal (samples)
        L,
        # start frequency (frequency at 0, normalized)
        f0,
        # end frequency (frequency reached at L, normalized)
        f1,
        # minimum
        x_min=0,
        # maximum
        x_max=1,
        # initial phase (between 0 and 1)
        phase=0.):
        # use like this:
        # chirp(...args...).sinusoid()
        self.n=np.arange(L)
        self.x=signal.chirp(self.n,f0,L,f1,phi=phase*360.)
        self.x_min=x_min
        self.x_max=x_max

    def sinusoid(self):
        x=self.x.copy()
        x+=1.
        x*=0.5
        x*=(self.x_max-self.x_min)
        x+=self.x_min
        return x

    def squarewave(self):
        y=np.zeros_like(self.x)
        y[self.x>=0]=self.x_max
        y[self.x<0]=self.x_min
        return y
