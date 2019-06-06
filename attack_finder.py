import librosa
import numpy as np
from scipy import signal

def find_attacks(
x,
# samples with values below this are not considered
thresh=1e-4,
# amount of smoothing for local maximum algorithm
smooth=0.9999,
# amount of samples to way before letting another attack through (limits the
# number of attacks that can occur over an interval of time)
N=1000):
    # find local maximum using full-wave rectified signal
    xm=fast_max(np.abs(x),alph=smooth)
    xmp=thresh_local_max_samples(np.diff(xm),alph=smooth,beta=1,N=N,thresh=thresh)
    return np.where(xmp > 0)[0]

def _shifted_lp(N=8,rs=80,Wn=0.5):
    b,a=signal.cheby2(N,rs,Wn)
    b=b*np.exp(complex("j")*Wn*np.pi*np.arange(len(b)))
    a=a*np.exp(complex("j")*Wn*np.pi*np.arange(len(a)))
    return (b,a)

def analytic_signal(x):
    """
    Make a rough analytic signal by filtering.
    """
    b,a=_shifted_lp(N=16)
    y=signal.lfilter(b,a,x)
    return y

def local_max(x,N=100):
    X=librosa.util.frame(x,frame_length=N,hop_length=1)
    y=np.max(X,axis=0)
    return y

def fast_max(x,alph=0.99):
    cmax=0
    ret=np.zeros_like(x)
    for n,x_ in enumerate(x):
        if x_ > cmax:
            cmax=x_
        ret[n]=cmax
        cmax*=alph
    return ret

def impulse_rate_limiter(x,alph=0.99,beta=0.5):
    cmax=0
    ret=np.zeros_like(x)
    for n,x_ in enumerate(x):
        if beta*x_ > cmax:
            cmax=x_
            ret[n]=x_
        cmax*=alph
    # first one is garbage
    ret[:2]=0
    return ret

def impulse_rate_limiter_samples(x,N=16000*.125):
    alph=0.999
    beta=np.power(alph,N)
    x=np.concatenate((np.ones((int(N)-1,)),x))
    ret=impulse_rate_limiter(x,alph=alph,beta=beta)
    return ret[int(N)-1:]

def thresh_local_max(x,alph=0.99,beta=0.5,alph2=.99):
    """
    This could be called "filter_attacks". What it does is accept plaubile
    attacks as impulses (perhaps from the diff of a local max or fast_max
    signal) and then output hopefully one for a group of attacks.
    """
    lm = fast_max(x,alph=alph)
    #lm=signal.lfilter([1],[1,-alph],x)
    ret=np.zeros_like(x)
    ret[x>=(lm*beta)]=1
    ret[0]=0
    #return ret
    return impulse_rate_limiter(ret,alph=alph2)
        
def thresh_local_max_samples(x,alph=0.99,beta=0.5,N=16000*.125,thresh=1e-3):
    """
    This could be called "filter_attacks". What it does is accept plaubile
    attacks as impulses (perhaps from the diff of a local max or fast_max
    signal) and then output hopefully one for a group of attacks.
    alph is how fast fast_max forgets the previous max
    beta is how much the max signal is scaled, impulses above this signal are
    not filtered
    N is the number of samples that samples are ignored from an initial sample
    signal values less than thresh are considered to be 0

    NOTE: This algorithm might anticipate attacks slightly too early
    """
    lm = fast_max(x,alph=alph)
    #lm=signal.lfilter([1],[1,-alph],x)
    ret=np.zeros_like(x)
    ret[(x>=(lm*beta))&(x>thresh)]=1
    ret[0]=0
    #return ret
    return impulse_rate_limiter_samples(ret,N=N)
