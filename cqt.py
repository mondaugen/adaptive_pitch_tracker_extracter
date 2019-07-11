import numpy as np

class cqt:
    """
    A class to carry out constant Q transforms.
    """
    def __init__(self,
        # A function that accepts a window length N and gives a window suitable
        # for Fourier analysis. E.g., you can use
        # get_window=lambda N: signal.get_window('hann',N,fftbins=True)
        get_window,
        # The length of the windows
        N,
        # The angular velocities (frequencies but in radians/s) at which to
        # compute exponentials which are then windowed. E.g., if you want to
        # project onto exponentials tuned to the notes from middle C to the C an
        # octave above, and the sampling rate is 48000, then pass
        # ws=[2*np.pi*440*np.power(2,(n-69)/12)/48000 for n in range(60,73)]
        ws):
        self.w=get_window(N)
        self.A=np.ndarray(shape=(len(ws),N),dtype='complex128')
        r=np.arange(N)
        for row,f in enumerate(ws):
            self.A[row,:]=np.exp(complex('j')*f*r)*self.w
    def __call__(self,x):
        return self.A @ x
            
