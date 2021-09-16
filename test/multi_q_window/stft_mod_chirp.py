import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import freq_dom_window

N_dft=2048
H_stft=512
window_type='blackman'
# Find the required Q of the lowest frequency
lowest_q=freq_dom_window.calc_Q_min(
1./N_dft,
freq_dom_window.window_types[window_type]['lobe_radius'],
N_dft)
N_v_q=5
v_q=[(v,q) for v,q in zip(np.linspace(0,0.5,N_v_q),np.linspace(1,4,N_v_q))]
#v_q=[(v,4) for v,q in zip(np.linspace(0,0.5,N_v_q),np.linspace(1,2,N_v_q))]
mqw=freq_dom_window.multi_q_window(N_dft,v_q,window_type)
# bin numbers
b=np.arange(1,N_dft//2)
# normalized bin frequencies
v=b/N_dft
# fourier transform transformation to apply special Q
R=mqw.R(v)
# length of signal
N_x=1000000
n_x=np.arange(N_x)
# starting chirp frequency
v0=0.25e-3
# ending chirp frequency
v1=0.5
# a chirp
x=signal.chirp(n_x,v0,N_x,v1)
v,t,Zxx=signal.stft(x,nperseg=N_dft,window='blackman',noverlap=N_dft-H_stft,return_onesided=False)
ZxxR=R@Zxx
N_show=100
for ws in mqw.Ws:
    plt.plot(np.arange(N_show),ws[:N_show])
plt.figure()
R_full=R.toarray()
for r in R_full:
    plt.plot(np.arange(len(r)),np.abs(r))
plt.matshow(20*np.log10(np.abs(ZxxR)))
#plt.matshow(20*np.log10(np.abs(Zxx)))

plt.show()
