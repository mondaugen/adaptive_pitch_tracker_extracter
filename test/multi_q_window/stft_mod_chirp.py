import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import freq_dom_window

N_dft=512
H_stft=128
mqw=freq_dom_window.multi_q_window(N_dft,[(0.005,1.15),(0.1,1.15)],'blackman')
# bin numbers
b=np.arange(N_dft)
# normalized bin frequencies
v=b/N_dft
# fourier transform transformation to apply special Q
R=mqw.R(v)
# length of signal
N_x=100000
n_x=np.arange(N_x)
# starting chirp frequency
v0=0.25e-3
# ending chirp frequency
v1=0.25
# a chirp
x=signal.chirp(n_x,v0,N_x,v1)
v,t,Zxx=signal.stft(x,nperseg=N_dft,window='blackman',noverlap=N_dft-H_stft)
ZxxR=R@Zxx
plt.matshow(20*np.log10(np.abs(ZxxR)))

plt.show()
