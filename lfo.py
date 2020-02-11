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
        self.phase=phase
        self.L=L
        self.f0=f0
        self.f1=f1
        self.n=np.arange(L)
        self.x=signal.chirp(self.n,self.f0,self.L,self.f1,phi=self.phase*360.)
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

    def sawtooth(self,falling=False):
        """
        if falling is False, ramps from self.x_min to self.x_max over one
        period, otherwise ramps from self.x_max to self.x_min.
        """
        x2=signal.chirp(self.n,self.f0,self.L,self.f1,
        phi=(self.phase*360.-90))
        s=np.arccos(self.x)/np.pi
        r=0.5*(1+(1-s)*np.power(-1,x2<0))
        if falling:
            return (self.x_max - self.x_min) * (1-r) + self.x_min
        else:
            return (self.x_max - self.x_min) * r + self.x_min

