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
R=interp1d([0,1],Rp,axis=0)(tR)

x=np.random.standard_normal(N)

iirlf=iirlf_f32(Rp.shape[1])
iirlfp=iirlf_f32_proc(x,R,1)
iirlf.proc(iirlfp)
_,_,_,Rest=ngal(iirlfp._out,3,alpha=0.001,beta=0.001)

for k,K in enumerate(Rest[1:].T):
    plt.plot(n,K,label="R_%d"%(k,))

plt.show()
