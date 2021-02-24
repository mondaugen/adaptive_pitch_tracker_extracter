import numpy as np
from scipy import signal
from common import frame
from exp_approx_line import exp_approx_line_coeffs

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
        # Buffer to hold past bin frequencies
        self.v_past=np.zeros((N,N_v))
        self.A0=np.zeros(self.N_v,dtype='complex128')
        self.Ap_plus=np.zeros((self.P-1,self.N_v),dtype='complex128')
        self.Ap_minus=np.zeros((self.P-1,self.N_v),dtype='complex128')
    def compute_dft(self,x,v):
        """
        Compute the next DFT centred on bins at normalized frequencies v.
        Generally v should be the same or close for every call for the result to
        be correct or close to correct.
        x is the value that will be shifted into self.buf (not in this function,
        but by an eventual call to self.shift_in)
        v must be a vector of length self.N_v.
        Returns Xv of length self.N_v which is the N-point DFT windowed by the
        sum-of-cosine window described by self.wp, centred on each bin in v.
        """
        # There are 2*self.P update equations if every cosine is split into 2
        # complex exponentials, however both of these exponentials are
        # everywhere 1 for p=0, so the there is one update equation for P=0 and
        # 2*(self.P-1) for the remaining self.P-1 coefficients
        self.A0=(np.exp(j*2*np.pi*v)*self.A0 
                + self.wp[0]*(np.exp(-j*2*np.pi*v*(self.N-1))*x 
                    - np.exp(j*2*np.pi*self.v_past[-self.N])*self.buf[-self.N]))
        for Ap, coef in zip([self.Ap_plus,self.Ap_minus],[-1,1]):
            for ap_k, (ap, wp) in enumerate(zip(Ap,self.wp[1:])):
                p=ap_k+1
                p_v=(p/self.N+coef*v)
                Ap[ap_k] = (np.exp(coef*j*2*np.pi*p_v)*ap
                    + wp*(np.exp(-coef*j*2*np.pi*p_v*(self.N-1))*x 
                        - np.exp(coef*j*2*np.pi*(p/self.N+coef*self.v_past[-self.N])
                            )*self.buf[-self.N]))
        # compute current DFT value
        Xv=self.A0+0.5*np.sum(self.Ap_plus+self.Ap_minus,axis=0)
        return Xv
    def shift_in(self,x):
        # shift in the current x and shift out self.buf[-N]
        # in faster implementation this would be done using a delay line or ring
        # buffer
        self.buf[1:]=self.buf[:-1]
        self.buf[-1]=x
    def shift_in_v(self,v):
        self.v_past[1:]=self.v_past[:-1]
        self.v_past[-1]=v
    def update(self,x,v):
        Xv=self.compute_dft(x,v)
        self.shift_in(x)
        self.shift_in_v(v)
        return Xv

class rec_dv_dft_sumcos(rec_dft_sumcos):
    """ Like rec_dft_sumcos but also returns an estimate of the derivative d/dv Xv. """
    def __init__(self,N,wp,N_v,max_err=1e-4):
        """
        max_err is (approximately) the maximum error tolerated when
        computing the ramp-approximating exponential.
        """
        rec_dft_sumcos.__init__(self,N,wp,N_v)
        self.B0_hat=np.zeros(self.N_v,dtype='complex128')
        self.Bp_hat_plus=np.zeros((self.P-1,self.N_v),dtype='complex128')
        self.Bp_hat_minus=np.zeros((self.P-1,self.N_v),dtype='complex128')
        alpha,beta,gamma,s0,s1=exp_approx_line_coeffs(self.N,self.N,max_err)
        self.alpha=alpha
        self.beta=beta
        self.gamma=gamma
        self.s0=s0
        self.s1=s1
    def compute_dv_dft(self,x,v):
        """
        Assumes self.A0, self.Ap_plus, and self.Ap_minus contain the computed
        terms from the most recent DFT, that is, assumes that compute_dft(x,v)
        has been called (with the same x and v) before this is called.
        Returns dvXv (of length self.N_v) which is the estimate of the
        derivative w.r.t. v of Xv at each bin in v.
        """
        self.B0_hat = (np.exp(-(self.alpha - j*2*np.pi*v))*self.B0_hat
            + self.wp[0]*np.exp(self.beta)*(
                np.exp((self.N-1)*(self.alpha-j*2*np.pi*v))*x
                - np.exp(-(self.alpha-j*2*np.pi*self.v_past[-self.N]))*self.buf[-self.N]))
        for Bp, coef in zip([self.Bp_hat_plus,self.Bp_hat_minus],[1,-1]):
            for bp_k, (bp, wp) in enumerate(zip(Bp,self.wp[1:])):
                p=bp_k+1
                p_v=(coef*p/self.N-v)
                Bp[bp_k] = (np.exp(-(self.alpha + j*2*np.pi*p_v))*bp
                    + wp*np.exp(self.beta)*(
                        np.exp((self.N-1)*(self.alpha+j*2*np.pi*p_v))*x
                        - np.exp(-(self.alpha+j*2*np.pi*(coef*p/self.N
                            -self.v_past[-self.N])))*self.buf[-self.N]))
        dvXv=(self.B0_hat + self.gamma*self.A0) + 0.5*np.sum(
            self.Bp_hat_plus + self.gamma * self.Ap_plus
            + self.Bp_hat_minus + self.gamma * self.Ap_minus,axis=0)
        return -j*2*np.pi*dvXv
    def update(self,x,v):
        Xv=self.compute_dft(x,v)
        dvXv=self.compute_dv_dft(x,v)
        self.shift_in(x)
        self.shift_in_v(v)
        return (Xv,dvXv)
