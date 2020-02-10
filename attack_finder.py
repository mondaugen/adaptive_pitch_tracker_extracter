import librosa
import numpy as np
from scipy import signal
import common
import spectral_difference

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
    xmp=thresh_local_max_samples(np.diff(xm),alph=smooth,beta=1,N=N,
        thresh=thresh)
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

def spectral_flux(x,H,W,window_type='hann'):
    """
    Frame-up x and take the DFT of each windowed frame. Compute the power
    spectrum for each frame. Form the sum of the absolute values of the
    differences between the adjacent bins. This is the spectral flux.
    Intuitively, a more jagged spectrum will give a greater value of spectral
    flux.
    """
    w=signal.get_window(window_type,W)
    w/=np.sum(W)
    spec_flux_i=0
    X=np.fft.fft(common.frame(x,H,W)*w[:,None],axis=0)
    spec_flux=np.sum(np.abs(np.diff(np.abs(X),axis=0)),axis=0)
    return spec_flux

def attack_region_estimation(x,H,W,alpha,thresh,window_type='boxcar'):
    """
    Return a pair (t,a) which has the times t and a signal a which is 1 if it is
    believed to be a part of an attack portion of a signal and 0 otherwise.
    H is the hop size, W the window size
    alpha is the amount of smoothing applied to the spectral flux before taking
    differences and thresholding to get the attack portions. 1 uses 100% of the
    current value and none of the past values 0.5 uses 50% of the current value
    and 50% of the last values, etc.
    window_type is the type of window used before taking the DFT.
    """
    spec_flux=spectral_flux(x,H,W,window_type)
    t_spec_flux=np.arange(len(spec_flux))*H+W/2
    smooth_filter_coeffs=([alpha],[1,-(1-alpha)])
    spec_flux=signal.lfilter(*smooth_filter_coeffs,spec_flux)
    spec_flux_diff=spec_flux[1:]-spec_flux[:-1]
    attacks=np.array(spec_flux_diff>thresh,dtype='float')
    return (t_spec_flux[1:],attacks)

def attacks_from_spectral_diff(
    # the signal whose attacks to estimate
    x,
    # hop size
    H=256,
    # window size
    W=1024,
    # window type 
    window_type='hann',
    # smoothing factor s. 0 < s <= 1
    # 1 means no smoothing, 0 would mean totally rejecting new values
    smoothing=1,
    # max filter discount rate, the time in samples until it reaches 1% of the
    # maximum
    lmax_filt_rate=16000,
    # the threshold in dB for the noise gate
    ng_th=-60):
    """ Returns pairs that are an estimation of the attack start and end times
    """

    # calculate the max filter rate
    # number of hops in the LMAX_FILT_RATE
    lmfr_n_H=lmax_filt_rate/H
    if lmfr_n_H <= 0:
        lmfr_n_H=1
    # what number to this power is .01 ?
    lmfr=np.power(.01,1/lmfr_n_H)

    sd=spectral_difference.spectral_diff(x,H,W,window_type)
    sd=spectral_difference.iir_avg(sd,smoothing)
    sd_maxs,sd_max_thresh=spectral_difference.discount_local_max(sd,lmfr)
    # set one_sided_max to 'left' so that the first point to become non-zero is
    # included as a local minimum
    sd_mins=spectral_difference.local_max(-sd,one_sided_max='left')
    sd_mins_filtered=spectral_difference.closest_index_after(sd_maxs,
        sd_mins,reverse=True)

    x_rms=spectral_difference.local_rms(x,H,W)
    sd_gate=(x_rms>np.power(10,ng_th/20)).astype(x.dtype)
    #sd_gate_t=np.arange(len(sd_gate))*H_SF/SAMPLE_RATE
    sd_maxs=spectral_difference.index_mask(sd_maxs,sd_gate)

    # add one hop to compensate for the differencing operation
    sd_maxs+=1
    sd_mins_filtered+=1
    
    # multiply by hop size because we want the indices in samples not hops
    sd_maxs*=H
    sd_mins_filtered*=H

    # compensate for window by offseting the maxima by half a window
    sd_maxs += W//2
    sd_mins_filtered += W//2

    return spectral_difference.up_down_match(sd_mins_filtered,sd_maxs)
