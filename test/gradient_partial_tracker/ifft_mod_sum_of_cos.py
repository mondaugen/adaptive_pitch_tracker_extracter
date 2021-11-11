from some_ft import mod_sum_of_cos_dft_dk, multi_mod_sum_of_cos_dft_dk
import numpy as np
import matplotlib.pyplot as plt

N=2048
L=(N)
W=L-1 # kinda arbitrary, smaller just gives oversampling
A=[0.42,0.5,0.08]
k=np.arange(N)
k0=10
X=mod_sum_of_cos_dft_dk(k0,k,A,L,W,N)
x=np.fft.ifft(X.astype('complex128'))
print(np.sum(np.abs(x.imag)))
plt.plot(k,x.real,label='real')
plt.plot(k,x.imag,label='imaginary')
plt.legend()
# partial multipliers
P=int(np.floor(0.5/(k0/N)))
p=(1+np.arange(P))
# amplitude scalars
B=1./p
dX_dk=multi_mod_sum_of_cos_dft_dk(B,k0*p,k,A,L,W,N)
dx_dk=np.fft.ifft(dX_dk.astype('complex128'))
plt.figure()
plt.plot(k,dx_dk.real,label='real')
plt.plot(k,dx_dk.imag,label='imaginary')
plt.legend()
plt.show()
