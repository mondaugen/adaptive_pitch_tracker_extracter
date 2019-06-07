# Find the attack of a signal
# At some point after this attack, make an estimate as to what decaying
# sinusoids are present in the signal
# Filter the signal with the inverse filters that would effectively remove these
# sinusoids. This gives the "error signal"
# We are interested in this signal because, together with the decaying sinusoids
# modelled as filters, we should be able to reconstruct the original signal
# Furthermore, we can change the centre frequencies of these filters, make their
# decays longer, etc.

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import librosa
import dft_freq_est
import attack_finder
from common import normalize

def calc_damped_sine_filters(bet,w,amps=None):
    """
    from amplitude decay coefficient bet and frequency w give the filter
    transfer function yielding a cosine with frequency w multiplied by a
    decaying exponential with decay bet multiplied by the unit step
    amps are the initial amplitudes of the signals
    """
    n_rows=len(bet)
    b=np.ones((n_rows,2))
    a=np.ones((n_rows,3))
    b[:,1]=-1*bet*np.cos(w)
    a[:,1]=-2*bet*np.cos(w)
    a[:,2]=bet*bet
    if amps is not None:
        #b=b*amps[:,None]
        g=amps*(1-bet*(1-bet))/(1-bet/2)
        b=b*g[:,None]
    return (b,a)

sr=16000
x=np.fromfile('/tmp/snd.f64')

N_W=1024
N_FFT=2048
N_H=512

peak_thresh=-100

# estimate where the attacks are
attack_i=attack_finder.find_attacks(x)
# how much time after the attack to wait before analysing
ad=int(0.25*sr)
ad+=attack_i[0]
# the longest decay possible (to account for analyses giving positive amplitude
# coefficients)
max_dec_s=10
max_bet=np.power(1e-3,1/(max_dec_s*sr))
print(max_bet)

# intialize the sinusoidal analyser
pa=dft_freq_est.peak_analyzer(
N_FFT,
N_W,
N_H)

# analyse some time after the attack
x_=x[ad:ad+pa.N_W+pa.H]
ph,amps,w_,bet=pa.freqs_amps(x_,max_n_peaks=50,peak_thresh=peak_thresh)
bet[bet>=1]=max_bet

pitch_trans=0
pitch_scale=np.power(2,pitch_trans/12.)

filt_b,filt_a=calc_damped_sine_filters(bet,w_)
filt_b_alt,filt_a_alt=calc_damped_sine_filters(bet,w_*pitch_scale)
err=x
for b,a in zip(filt_b,filt_a):
    err=signal.lfilter(a,b,err)

y=err
for b,a in zip(filt_b_alt,filt_a_alt):
    y=signal.lfilter(b,a,y)
y=normalize(y)

imp=np.zeros_like(x)
imp[0]=1
for b,a in zip(filt_b_alt,filt_a_alt):
    imp=signal.lfilter(b,a,imp)
imp=normalize(imp)

err_re=signal.resample(err,int(len(err)/pitch_scale))
for b,a in zip(filt_b_alt,filt_a_alt):
    err_re=signal.lfilter(b,a,err_re)
err_re=normalize(err_re)

err.tofile('/tmp/err.f64')
y.tofile('/tmp/y.f64')
imp.tofile('/tmp/imp.f64')
err_re.tofile('/tmp/err_re.f64')
