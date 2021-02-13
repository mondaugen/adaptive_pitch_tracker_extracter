import numpy as np
from scipy import signal
from common import frame
# Test the recursive implementation of a hann windowed sequence
j=complex('j')
class rec_dft_sumcos:
    def __init__(self,N,wp,N_v):
        """
        Initialize the recursive DFT of a sequence windowed by a sum-of-cosines window.
        N is the length of the window.
        wp is the p window coefficients, e.g., for a hann window these are
        [0.5,-0.5] (notice that the window starts at n = 0 and is centred on
        N/2).
        N_v is the number of bins to transform.
        """
        self.N=N
        self.wp=wp
        self.P=len(self.wp)
        # Buffer to hold past signal values, in a faster implementation this
        # would be a ring buffer or a delay line.
        self.buf=np.zeros(N)
        self.N_v=N_v
        self.A0=np.zeros(self.N_v,dtype='complex128')
        self.Ap_plus=np.zeros((self.P-1,self.N_v),dtype='complex128')
        self.Ap_minus=np.zeros((self.P-1,self.N_v),dtype='complex128')
    def update(self,x,v):
        """
        Compute the next DFT centred on bin v. Generally v should be the same or
        close for every call for the result to be correct or close to correct.
        v must be a vector of length self.N_v.
        """
        # There are 2*self.P update equations if every cosine is split into 2
        # complex exponentials, however both of these exponentials are
        # everywhere 1 for p=0, so the there is one update equation for P=0 and
        # 2*(self.P-1) for the remaining self.P-1 coefficients
        self.A0=(np.exp(j*2*np.pi*v)*self.A0 
                + self.wp[0]*(np.exp(-j*2*np.pi*v*(self.N-1))*x 
                    - np.exp(j*2*np.pi*v)*self.buf[-self.N]))
        for Ap, coef in zip([self.Ap_plus,self.Ap_minus],[-1,1]):
            for ap_k, (ap, wp) in enumerate(zip(Ap,self.wp[1:])):
                p=ap_k+1
                p_v=(p/self.N+coef*v)
                Ap[ap_k] = (np.exp(coef*j*2*np.pi*p_v)*ap
                    + 0.5*wp*(np.exp(-coef*j*2*np.pi*p_v*(self.N-1))*x 
                        - np.exp(coef*j*2*np.pi*p_v)*self.buf[-self.N]))
        # shift in the current x and shift out self.buf[-N]
        self.buf[1:]=self.buf[:-1]
        self.buf[-1]=x
        # compute current DFT value
        Xv=self.A0+np.sum(self.Ap_plus+self.Ap_minus,axis=0)
        return Xv

if __name__ == '__main__':
    N_x=10000
    N_w=1024
    w=signal.get_window('hann',N_w)
    x=np.random.standard_normal(N_x)
    x_slow=np.concatenate((np.zeros(N_w-1),x))
    wp=np.array([0.5,-0.5])
    v=np.array([0.001,0.002])
    rec_def=rec_dft_sumcos(N_w,wp,len(v))
    rec_def_slow=np.exp(-j*2*np.pi*v[:,None]*np.arange(N_w))*w
    print('rec_def_slow',rec_def_slow)
    x_framed=frame(x_slow,1,N_w)
    print('x_framed[:,0]',x_framed[:,0])
    print('x[0]',x[0])
    print('x_framed.shape',x_framed.shape)
    Xv_slow=rec_def_slow@x_framed
    Xv=np.concatenate([rec_def.update(x_,v)[:,None] for x_ in x],axis=1)
    print('Xv_slow.shape',Xv_slow.shape)
    print('Xv.shape',Xv.shape)
    Xvs=Xv_slow[:,:10].T
    Xvf=Xv[:,:10].T
    print('Xv_slow:',Xvs)
    print('Xv:',Xvf)
    print('Xv_slow/Xv:',np.abs(Xvs/Xvf))
