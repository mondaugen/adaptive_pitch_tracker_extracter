import numpy as np
import matplotlib.pyplot as plt
from some_ft import dirichlet, dirichlet_dk
from some_sig import rectangular
from dftdk import dft_dk, dft_bin, ps_dk

print_ = print

show_plot=True
print_data=True

def print(*args):
    if print_data:
        print_(*args)

j=complex('j')
N=4096
n=np.arange(N)
W=7
k=np.linspace(0,100,123) # non-integer bins
x=rectangular(n,W,N)
plt.plot(n,x)
X=np.fft.fft(x)
Xni=np.real(dft_bin(x,k))
Xth=dirichlet(n,W,N)
plt.figure()
plt.plot(n,np.real(X),label='computed')
plt.plot(n,Xth,label='theoretical')
plt.plot(k,Xni,'.',label='non-integer k')
plt.legend()
plt.figure()
dXth=dirichlet_dk(n,W,N)
dXthfd=np.diff(Xth)
n_=np.concatenate((np.arange(N//2),np.arange(N//2)-N//2))
nx=n_*x
print("nx*W:")
print(nx*W)
dX=dft_dk(N)(x)#np.fft.fft(-j*2*np.pi/N*nx)
dXni=np.real(dft_bin(x,k,p=1))
plt.plot(n,np.real(dX),label='computed_real')
plt.plot(n,np.imag(dX),label='computed_imag')
plt.plot(n,np.abs(dX),'--',label='computed_abs')
plt.plot(n[:-1]+0.5,dXthfd,'--',label='finite difference')
plt.plot(n,dXth,'--',label='theoretical')
plt.plot(k,dXni,'.',label='non-integer k')
plt.legend()
print("sum(abs(computed_real-theoretical)):",
np.sum(np.abs(np.real(dX)-dXth)))

plt.figure()
# estimate of second derivative
dX2fd=np.diff(np.real(dX))
# computed with fourier transform
dX2=dft_dk(N,p=2)(x)
plt.title('2nd derivative')
plt.plot(n[:-1]+0.5,dX2fd,label='fininte difference')
plt.plot(n,np.real(dX2),label='computed_real')
plt.legend()

plt.figure()
# power spectrum
ps=np.real(X*np.conj(X))
plt.title('power spectrum')
plt.plot(n,ps)

plt.figure()
# derivative of power spectrum
dPfd=np.diff(ps)
dP=ps_dk(X,dX)
plt.title("power spectrum derivative")
plt.plot(n[:-1]+0.5,dPfd,label='finite difference')
plt.plot(n,dP,label='computed')
plt.plot(k,ps_dk(Xni,dXni,extract=np.real),label='non-integer k real')
plt.plot(k,ps_dk(Xni,dXni,extract=np.imag),label='non-integer k imaginary')
plt.legend()

# oversampled using zero-padding
x_zp = np.concatenate((x[:N//2],np.zeros(N),x[N//2:]))
X_zp = np.fft.fft(x_zp)
n_zp = np.arange(0,N,0.5)
plt.figure()
plt.plot(n,np.real(X),label="original")
plt.plot(n_zp,np.real(X_zp),label="zero-padded")
plt.legend()


if show_plot:
    plt.show()
