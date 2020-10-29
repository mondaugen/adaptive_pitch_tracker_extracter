import numpy as np
from scipy import signal, interpolate
import common
import subprocess

def swap(a,b):
    return b,a

def a_to_r(a):
    """ Convert feedback filter coefficients to reflection coefficients. This assumes you
    include the coefficient for y[n], where it divides all the coefficients by
    it (so the first coefficient is 1), then discards the y[n] coefficient
    (because it knows it is 1). """
    a=a.copy()
    a=a/a[0]
    a=a[1:]
    p=len(a)
    r=np.zeros_like(a)
    r[-1]=a[-1]
    for j in np.arange(p-1)[::-1]:
        s=1/(1-r[j+1]*np.conj(r[j+1]))
        a[:j+1]=s*(a[:j+1]-r[j+1]*np.conj(a[:j+1][::-1]))
        r[j]=a[j]
    return r

def r_to_A(r):
    """
    Convert reflection coefficients to feedback coefficients. This returns all
    the intermediate filter orders along the rows, so the first row will be
    [1,0,0, .., 0], the second row [1,a_1[1], 0, 0, ..., 0], etc. So if you just
    want the highest filter order obtainable with these r coefficients just take
    A[-1,:] of the resulting matrix A.
    """
    P=len(r)
    A=np.zeros((P+1,P+1),dtype=r.dtype)
    A[0,0]=1
    for i in np.arange(1,P+1):
        A[i,:i+1]=A[i-1,:i+1]+r[i-1]*np.conj(A[i-1,:i+1][::-1])
    return A

def b_A_to_c(b,A):
    """
    Convert feedforward filter coefficients to 'C coefficients'. See p. 306 in
    Hayes book. These are used for implementing a lattice filter with poles and
    zeros.
    """
    assert(len(b)==A.shape[0])
    # This is to match the notation in Hayes
    q=len(b)-1
    c=np.zeros_like(b)
    c[q]=b[q]
    for k in np.arange(q)[::-1]:
        c[k]=b[k]-np.sum(c[k+1:q+1]*np.conj(np.diag(A,-k)[1:]))
    return c

def c_A_to_b(c,A):
    """ Convert 'C coefficients' to feedforward coefficients. """
    b=np.zeros_like(c)
    q=len(b)-1
    for k in np.arange(q+1):
        b[k]=np.sum(c[k:q+1]*np.conj(np.diag(A,-k)))
    return b

def c_r_to_b_a(c,r):
    A=r_to_A(r)
    b=c_A_to_b(c,A)
    return (b,A[-1,:])

def b_a_to_r_c(b,a):
    r=a_to_r(a)
    A=r_to_A(r)
    c=b_A_to_c(b,A)
    return (r,c)

# burg and burg2 give good estimate of coefficients, but I don't know how to
# interpret the error signal because it does not match the original signal
# filtered by the FIR given by the coefficients
def burg(x,p):
    """
    The Burg algorithm
    returns the p reflection coefficients and the pth forward prediction error
    (the error signal), estimated from signal x using Burg's method.
    """
    N=len(x)
    temp1=x
    temp2=x
    gamma=np.zeros(p,dtype=x.dtype)
    for j in range(p):
        eplus = temp1[1:]
        eminus = temp2[:-1]
        gamma[j] = (-2*np.sum(np.conj(eminus)*eplus)
            / (np.sum(np.conj(eplus)*eplus)+np.sum(np.conj(eminus)*eminus)))
        temp1 = eplus + gamma[j] * eminus
        temp2 = eminus + np.conj(gamma[j])*eplus
    return gamma, eplus, eminus

def burg2(x,p):
    N=len(x)
    x=np.concatenate((np.zeros(p),x))
    temp1=x
    temp2=x
    gamma=np.zeros(p,dtype=x.dtype)
    for j in range(p):
        eplus = temp1[1:]
        eminus = temp2[:-1]
        gamma[j] = (-2*np.sum(np.conj(eminus)*eplus)
            / (np.sum(np.conj(eplus)*eplus)+np.sum(np.conj(eminus)*eminus)))
        temp1 = eplus + gamma[j] * eminus
        temp2 = eminus + np.conj(gamma[j])*eplus
    return gamma, eplus, eminus

# This gives an interpretable error signal but the coefficient estimates are very poor
def burgN(x,p):
    N=len(x)
    ej_1p=np.zeros(N)
    ejp  =np.zeros(N)
    ej_1m=np.zeros(N)
    ejm=np.zeros(N)
    G=np.zeros(p)
    E=np.zeros(p)
    ej_1p[:]=x
    ej_1m[1:]=x[:-1]
    for j in range(p):
        G[j] = 2*np.sum(ej_1p*np.conj(ej_1m)) / np.sum(np.abs(ej_1p)**2+np.abs(ej_1m)**2)
        ejp[:] = ej_1p + G[j] * ej_1m
        ejm[1:] = ej_1m[:-1] + np.conj(G[j]) * ej_1p[:-1]
        ejmN = ej_1m[-1] + np.conj(G[j]) * ej_1p[-1]
        ejp,ej_1p = swap(ejp,ej_1p)
        ejm,ej_1m = swap(ejm,ej_1m)
    return G,ej_1p,ej_1m,E



def lattice_filter_proc(x,R,C):
    """
    R and C have the columns filled with the coefficients for a specific time
    step and time advances along the rows.
    """
    # This works becase len always gets length of 1st dimension
    P = len(R)
    in_file_name=common.mktemp(suf='.f32')
    out_file_name=common.mktemp(suf='.f32')
    R_file_name=common.mktemp(suf='.f32')
    C_file_name=common.mktemp(suf='.f32')
    x.astype('float32').tofile(in_file_name)
    R.T.astype('float32').tofile(R_file_name)
    C.T.astype('float32').tofile(C_file_name)
    env=dict(
        P='%d'%(P,),
        IN_PATH=in_file_name,
        OUT_PATH=out_file_name,
        R_PATH=R_file_name,
        C_PATH=C_file_name)
    subprocess.run('src/test/bin/lattice_filter_proc',env=env)
    y=np.fromfile(out_file_name,dtype='float32')
    return y

def iir_lattice_filter_proc(x,R):
    # TODO This doesn't work, let's make one that does
    P = len(R)
    b = np.zeros
    b[0] = 1
    c = b_A_to_c(b)
    C = np.zeros((P+1,len(x)))
    C[:P,:]=c
    return lattice_filter_proc(x,R,C)

class filter_interp_table:

    def __init__(self,
        # function or callable object accepting filter order n and critical
        # frequency f returning a tuple (b,a) like the scipy.signal filter
        # design functions.
        # frequency is normalized between 0 and 1 where 1 is the sampling rate
        design_fun,
        # filter_order+1 is the number of numerator and denominator coefficients
        filter_order,
        # minimum frequency: filters are designed at exponentially spaced
        # frequencies min_freq, min_freq*2**(1/designs_per_octave),
        # min_freq*2**(2/designs_per_octave), ...
        # minimum frequency is normalized between 0 and 1 where 1 is sampling rate
        min_f,
        # number of designs per octave: the spacing between known filter designs
        # (the designs for other frequencies are obtained by interpolating the
        # coefficients in a stable way)
        designs_per_octave=1,
        # maximum frequency: filter designs stop at this frequency
        max_f=0.5,
        dtype=np.float64,
        # If none, the out_of_bounds_values are just fixed at the left or right
        # endpoint of the known values
        # if 'lowpass' the filter response corresponding to a gain of 0 is put
        # at frequency 0 and filter response of gain of 1 for frequency = max_f
        # if 'highpass' the opposite is done (gain of 1 for frequency 0, gain of
        # 0 for max_f)
        out_of_bounds_values=None,
        check_filter_stability=True
        ):
        """
        TIPS:
        It is a good idea to set check_filter_stability=True.
        If you are getting unstable filters, try higher and lower
        values of min_f and max_f respectively. Ultimately this is the fault of
        the filter design routine though haha.
        NOTE:
        It has been observed with out_of_bounds_values='lowpass' and a min_f of
        0.003 (around 50Hz at 16KHz sample rate) that a subtle pop (not a click)
        occurs at the end of a ramp to a frequency of 0. This seems to be
        because the reflection coefficients were far from 0 at min_f (in fact
        they were very nearly 1 or -1), but to give a gain of zero, all but one
        has to ramp to 0. This slow ramping seems to cause a ringing in the
        filter. Setting a lower min_f made the ringing oscillate at an inaudible
        frequency, but this might not be possible for higher order designs. A
        way to get a lower min_f is to have the filter design function return
        smaller filter orders closer to min_f. Finally you could just ramp to a
        minimum of min_f (set out_of_bounds_values=None) and apply a little
        amplitude ramp.
        """
        assert(filter_order>0)
        assert(min_f>0)
        assert(max_f>0)
        assert(designs_per_octave>0)
        n_lim=int(np.floor(designs_per_octave*np.log2(max_f/min_f)))
        self.R=np.zeros((filter_order,n_lim),dtype=dtype)
        self.C=np.zeros((filter_order+1,n_lim),dtype=dtype)
        self.design_f=min_f*np.power(2,np.arange(n_lim)/designs_per_octave)
        self.max_f=max_f
        self.min_f=min_f
        self.filter_order=filter_order
        self.dtype=dtype
        for n,f in enumerate(self.design_f):
            b,a=design_fun(filter_order,f)
            r,c=b_a_to_r_c(b,a)
            if check_filter_stability:
                if np.any(np.abs(r)>1):
                    raise ValueError('Filter unstable for design at frequency %f' % (f,))
            self.R[:,n]=r
            self.C[:,n]=c

        if out_of_bounds_values is not None:
            self.R=np.concatenate((
                np.zeros((self.filter_order,1),dtype=self.dtype),
                self.R,
                np.zeros((self.filter_order,1),dtype=self.dtype)),axis=1)
            self.C=np.concatenate((
                np.zeros((self.filter_order+1,1),dtype=self.dtype),
                self.C,
                np.zeros((self.filter_order+1,1),dtype=self.dtype)),axis=1)
            self.R[0,0]=1
            self.R[0,-1]=1
            self.design_f=np.concatenate(([0],self.design_f,[self.max_f]))
            a=np.zeros(filter_order+1,dtype=dtype)
            a[0]=1
            b=np.zeros(filter_order+1,dtype=dtype)
            r_nulling,c_nulling=b_a_to_r_c(b,a)
            b[0]=1
            r_unity,c_unity=b_a_to_r_c(b,a)
            if out_of_bounds_values == 'lowpass':
                self.R[:,0]=r_nulling
                self.C[:,0]=c_nulling
                self.R[:,-1]=r_unity
                self.C[:,-1]=c_unity
            elif out_of_bounds_values == 'highpass':
                self.R[:,-1]=r_nulling
                self.C[:,-1]=c_nulling
                self.R[:,0]=r_unity
                self.C[:,0]=c_unity

        self.R_interp=interpolate.interp1d(
            self.design_f,
            self.R,
            fill_value=(self.R[:,0],self.R[:,-1]),
            bounds_error=False)
        self.C_interp=interpolate.interp1d(
            self.design_f,
            self.C,
            fill_value=(self.C[:,0],self.C[:,-1]),
            bounds_error=False)

    def lookup(self,f):
        # Returns reflection coefficients R and C coefficients for each value of
        # f interpolated from values in the table
        R=self.R_interp(f)
        C=self.C_interp(f)
        return (R,C)
