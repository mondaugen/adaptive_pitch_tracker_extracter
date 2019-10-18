# A function giving estimation of attacks using spectral difference

import numpy as np
from scipy import signal
import common
import matplotlib.pyplot as plt

def spectral_diff(x,H,W,window_type):
    w=signal.get_window(window_type,W)
    x_framed=common.frame(x,H,W)*w[:,None]
    X_framed=np.fft.rfft(x_framed,axis=0)/np.sum(W)
    xd_abs=np.abs(X_framed)
    sd=xd_abs[:,1:]-xd_abs[:,:-1]
    sd[sd<0]=0
    sd=np.sum(sd,axis=0)
    return sd

def high_frequency_weight(x,H,W,window_type):
    w=signal.get_window(window_type,W)
    x_framed=common.frame(x,H,W)*w[:,None]
    X_framed=np.fft.rfft(x_framed,axis=0)/np.sum(W)
    x_abs=np.abs(X_framed)
    scalar=np.arange(x_abs.shape[0])[:,None]
    hfw=np.sum(x_abs*scalar,axis=0)
    return hfw

def local_rms(x,H,W):
    x_framed=common.frame(x,H,W)
    x_rms=np.sqrt(np.sum(x_framed**2,axis=0))
    return x_rms

def local_max(x):
    """ for vectors only """
    gtr=np.concatenate((x[:-1]>=x[1:],np.zeros((1),dtype='bool')),axis=0)
    gtl=np.concatenate((np.zeros((1),dtype='bool'),x[1:]>x[:-1]),axis=0)
    return np.where(gtr&gtl)[0]

def local_max_mat(x):
    """ for matrices, finds the local maxima within the columns """
    gtr=np.concatenate((x[:-1,:]>=x[1:,:],np.zeros((1,x.shape[1]),dtype='bool')),axis=0)
    gtl=np.concatenate((np.zeros((1,x.shape[1]),dtype='bool'),x[1:,:]>x[:-1,:]),axis=0)
    return np.where(gtr&gtl)

def filtered_local_max(x,H,W,a):
    """
    a value of x is deemed the local maximum if it is a local maximum within
    a window of size W and is greater than a times the minimum local maximum
    """
    x_f=common.frame(x,H,W)
    res=np.zeros(len(x))
    for h,c in zip(np.arange(0,len(x)-W,H),x_f.T):
        max_n=local_max(c)
        if (len(max_n) == 0):
            # no local maxima
            continue
        if len(max_n) == 1:
            # there's no other local maximum to compare with, so we compare with
            # the mean value of the frame
            if c[max_n[0]] > a*np.mean(c):
                res[max_n[0]+h] = 1
            continue
        c_max=c[max_n]
        max_n_sorted=np.argsort(c_max)
        #if (c_max[max_n_sorted[-1]] > a*c_max[max_n_sorted[0]]):
        if (c_max[max_n_sorted[-1]] > a*np.mean(c)):
            res[max_n[max_n_sorted[-1]]+h] = 1
    return np.where(res>0)[0]

INPUT=common.get_env('INPUT',check_if_none=True)
WINDOW_TYPE=common.get_env('WINDOW_TYPE',default='hann')
SAMPLE_RATE=common.get_env('SAMPLE_RATE',conv=float,default=16e3)
H=common.get_env('H',conv=int,default=2)
W=common.get_env('W',conv=int,default=8)
H_RMS=common.get_env('H_RMS',conv=int,default=512)
W_RMS=common.get_env('W_RMS',conv=int,default=2048)
H_LMAX=common.get_env('H_LMAX',conv=int,default=5)
W_LMAX=common.get_env('W_LMAX',conv=int,default=10)
# How many times greater the max has to be to get counted
A_LMAX=common.get_env('A_LMAX',conv=float,default=1.5)

X_LIM=common.get_env('X_LIM',conv=eval,default=(0,2))

x=np.fromfile(INPUT)
t_x=np.arange(len(x))/SAMPLE_RATE
sd=spectral_diff(x,H,W,WINDOW_TYPE)
n_h=np.arange(H,len(x)-W,H)
t_sd=n_h[:-1]/SAMPLE_RATE
t_h=n_h/SAMPLE_RATE
hfw=high_frequency_weight(x,H,W,WINDOW_TYPE)
print(np.any(np.isnan(hfw)))

n_h_rms=np.arange(H_RMS,len(x)-W_RMS,H_RMS)
t_h_rms=n_h_rms/SAMPLE_RATE

x_rms=local_rms(x,H_RMS,W_RMS)

fig,axs=plt.subplots(4,1)
axs[0].plot(t_x,x)
axs[0].set_title('signal')
sd_log=np.log(sd+1e-6)
sd_l_max=filtered_local_max(sd,H_LMAX,W_LMAX,A_LMAX)
axs[1].plot(t_sd,sd_log)
axs[1].plot(t_sd[sd_l_max],sd_log[sd_l_max],'.')
axs[1].set_title('spectral difference')
hfw_l_max=filtered_local_max(hfw,H_LMAX,W_LMAX,A_LMAX)
axs[2].plot(t_h,hfw)
axs[2].plot(t_h[hfw_l_max],hfw[hfw_l_max],'.')
axs[2].set_title('high frequency weighting')
axs[3].plot(t_h_rms,x_rms)
axs[3].set_title('rms')

for ax in axs:
    ax.set_xlim(*X_LIM)

plt.tight_layout()
plt.show()

