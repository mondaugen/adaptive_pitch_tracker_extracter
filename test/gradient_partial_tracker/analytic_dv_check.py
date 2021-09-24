import numpy as np
import matplotlib.pyplot as plt

def dirichlet(k,W,N):
    ret = np.sin((np.pi*W*k)/N)/(W*np.sin(np.pi*k/N))
    ret[k == 0] = 1
    return ret

def dirichlet_dk(k,W,N):
    ret = np.pi/(N*np.sin(np.pi*k/N))*(np.cos(np.pi*W*k/N) - np.sin(np.pi*W*k/N)*np.cos(np.pi*k/N)/(W*np.sin(np.pi*k/N)))
    ret[k == 0] = 0
    return ret

def rectangular(n,W,N):
    ret = np.zeros(len(n))
    ret[(2*n < W)|(2*(N-n)<W)] = 1./W
    return ret
j=complex('j')
N=256
n=np.arange(N)
W=7
x=rectangular(n,W,N)
plt.plot(n,x)
X=np.fft.fft(x)
Xth=dirichlet(n,W,N)
plt.figure()
plt.plot(n,np.real(X),label='computed')
plt.plot(n,Xth,label='theoretical')
plt.legend()
plt.figure()
dXth=dirichlet_dk(n,W,N)
dXthfd=np.diff(Xth)
n_=np.concatenate((np.arange(N//2),np.arange(N//2)-N))
# This I got to by guessing but it's getting closer...
#dX=-j*2*np.pi/N*np.fft.fft(n_*x/N)
dX=np.fft.fft(-j*2*np.pi/N*n*x)
plt.plot(n,np.real(dX),label='computed_real')
plt.plot(n,np.imag(dX),label='computed_imag')
plt.plot(n,np.abs(dX),label='computed_abs')
plt.plot(n[:-1]+0.5,dXthfd,label='finite difference')
plt.plot(n,dXth,label='theoretical')
print("quotient")
print(np.angle(dX))
#print("abs quotient")
#print(np.abs(dX/dXth))
plt.legend()
plt.show()
