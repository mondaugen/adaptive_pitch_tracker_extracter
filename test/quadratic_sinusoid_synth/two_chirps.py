import numpy as np
import cubic_sinusoid_synth as cssyn
import matplotlib.pyplot as plt
from scipy import signal
from common import normalize


j=complex('j')

n_partials=5
theta=0.1*(np.arange(n_partials)+1)
F=100
omega=np.multiply.outer(
    np.concatenate((
        np.linspace(1,2,F),
        np.linspace(2,0.125,F),
        np.linspace(0.125,0.2,3*F),
        )),
    0.1*(np.arange(n_partials)+1)
)
H=256
phs=cssyn.quadratic_phase_poly_interp(H,theta,omega)
x=np.sum(np.exp(j*phs),axis=0)

print(x.shape)
N=2048
N_h=512
plt.specgram(x,NFFT=N,Fs=1,
window=signal.get_window('hann',N),
noverlap=N-N_h,sides='onesided')

normalize(np.real(x).astype('float64')).tofile('/tmp/x.f64')

plt.figure()

n=np.arange(len(x))
plt.plot(n,x)

plt.figure()

n=np.arange(phs.shape[1])
plt.plot(n,phs.T)

plt.show()

