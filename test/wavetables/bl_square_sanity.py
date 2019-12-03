# There shouldn't be ringing in a bandlimited square wave, no?
import numpy as np
import matplotlib.pyplot as plt

SR=1
N_dft=8192
f0_dft=SR/N_dft
a=f0_dft/0.01
print(a)
X=2*(np.sinc(np.arange(N_dft)*f0_dft/a)/a).astype('complex128')
X*=(np.arange(N_dft)*f0_dft<0.005).astype('float')
X[0]=0
x=np.real(np.fft.ifft(X)).astype('float32')

P=100
x.tofile('/tmp/fourier_synth.f32')
with open('/tmp/fourier_synth.f32','a') as f:
    for p in range(P):
        x.tofile(f)

plt.plot(np.arange(N_dft),x)
plt.show()

