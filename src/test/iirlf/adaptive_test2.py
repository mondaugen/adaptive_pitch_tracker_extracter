import numpy as np
from scipy.interpolate import interp1d
from iir_lattice_filter import iirlf_f32, iirlf_f32_proc
from gal_alexander import ngal
import matplotlib.pyplot as plt

# Filter a signal with time-varying reflection coefficients
# Try to extract the coefficients using an adaptive filter

N=100000
n=np.arange(N)

Rp=np.array([
    [0.75,0.25,-0.5],
    [-0.75,-0.25,0.5],
])

tR=np.linspace(0,1,N)
R=interp1d([0,1],Rp,axis=0,kind='linear')(tR)

x=np.random.standard_normal(N)

iirlf=iirlf_f32(Rp.shape[1])
iirlfp=iirlf_f32_proc(x,R,1)
iirlf.proc(iirlfp)
Ef,_,_,Rest=ngal(iirlfp._out,3,alpha=0.3,beta=0.001,normalize=False)

iirlf2=iirlf_f32(Rest.shape[1])
iirlfp2=iirlf_f32_proc(Ef[:,-1],Rest,1)
iirlf2.proc(iirlfp2)
yy=iirlfp2._out

fig,ax=plt.subplots(3,1)
for k,K in enumerate(Rest[1:].T):
    ax[0].plot(n,K,label="R_%d"%(k,))
ax[1].plot(n,Ef[1:,-1])
ax[2].plot(n,iirlfp._out,label='original')
ax[2].plot(n,yy[1:],label='reconstructed')
ax[2].legend()
print(np.mean(np.abs(iirlfp._out-yy[1:])))

plt.show()
