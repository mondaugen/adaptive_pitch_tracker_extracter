# Find the average derivative of the DFT w.r.t k (DFT(k)/dk) for harmonically
# spaced k on a signal that is itself a sum of harmonically spaced complex
# exponentials
# Then confirm that the same can be done using a signal that is simply the sum
# of the exponentials at the bin frequencies k (well -k because they are fourier
# transform kernels)

import numpy as np
from some_ft import (multi_mod_sum_of_cos_dft_k,
                     multi_mod_sum_of_cos_dft_dk,
                     multi_mod_sum_of_cos_dft_k,
                     normalize_sum_of_cos_A)
from some_sig import sum_of_cos, mod_sum_of_cos, multi_mod_sum_of_cos, nonsym_multi_mod_sum_of_cos
import dftdk
from dftdk import multiply_ramp


N=2048
W=1023
L=W+1
k0_sig=N*0.1
k0_anl=N*np.power(0.1,1/12) # can specify in semitones
v0_sig=k0_sig/N
v0_anl=k0_anl/N
# partial multipliers
P=10
p=(1+np.arange(P))
# amplitude scalars
B=1./p
raw_A=[0.42,0.5,0.08]
A=normalize_sum_of_cos_A(raw_A,L,W,N)
dX_dk=multi_mod_sum_of_cos_dft_dk(B,k0_sig*p,k0_anl*p,A,L,W,N)
avg_dX_dk=dX_dk.mean()

# The average derivative
print("avg_dX_dk:",avg_dX_dk)

# The value of the DFT at each analysis bin
X=multi_mod_sum_of_cos_dft_k(B,k0_sig*p,k0_anl*p,A,L,W,N)
X_sum=np.sum(X)
print("X:",X)
print("sum(X):",X_sum)

# Now the DFT done directly using two time-domain signals
x_anl=multi_mod_sum_of_cos(np.ones_like(B),v0_anl*p,A,L,W,N)
x_sig=multi_mod_sum_of_cos(B,v0_sig*p,[1.],L,W,N)
X_sum_td=np.sum(np.conj(x_anl)*x_sig)
print("X_sum_td:",X_sum_td)
print("X_sum_td/X_sum:",X_sum_td/X_sum)

# DFT dk using time-domain signal
sig=multi_mod_sum_of_cos(B,v0_sig*p,[1.],L,W,N)
#dX_dk_kern=multi_mod_sum_of_cos(np.ones_like(B)/P,v0_anl*p,A,L,W,N)
# Done by multiplying the harmonically spaced sinusoids with the window after
# summing them to prove it can be done this way
dX_dk_kern=multi_mod_sum_of_cos(
    np.ones_like(B)/P,v0_anl*p,[1],L,W,N)*sum_of_cos(A,L,W,N)
avg_dX_dk_td=(np.conj(dX_dk_kern)*multiply_ramp(sig,N,1)).sum()
print("avg_dX_dk_td:",avg_dX_dk_td)
print("avg_dX_dk_td/avg_dX_dk:",avg_dX_dk_td/avg_dX_dk)

# DFT dk using the class
hgtd=dftdk.harm_grad_td(raw_A,L,W,N,
    harm_sig=dftdk.harm_grad_td.default_harm_sig(B=np.ones_like(B)/P)
)
avg_dX_dk_td_hgtd=hgtd.dX_dk(sig,k0_anl)
print("avg_dX_dk_td_hgtd:",avg_dX_dk_td_hgtd)
print("avg_dX_dk_td_hgtd/avg_dX_dk:",avg_dX_dk_td_hgtd/avg_dX_dk)
avg_dX_dk_td_hgtd_fd=hgtd.dX_dk_fd(sig,k0_anl)
print("avg_dX_dk_td_hgtd_fd:",avg_dX_dk_td_hgtd_fd)
# non-symmetrical signal, therefore complex fourier transform
sig_nonsym=nonsym_multi_mod_sum_of_cos(B,v0_sig*p,N)
avg_d_log_ps_nonsym=hgtd.d_log_ps_dk(sig,k0_anl)
avg_d_log_ps_nonsym_fd=hgtd.d_log_ps_dk_fd(sig,k0_anl)
print("avg_d_log_ps_nonsym:",avg_d_log_ps_nonsym)
print("avg_d_log_ps_nonsym_fd:",avg_d_log_ps_nonsym_fd)
print("avg_d_log_ps_nonsym/avg_d_log_ps_nonsym_fd:",avg_d_log_ps_nonsym/avg_d_log_ps_nonsym_fd)
# TODO: It seems the finite differences are wrong !

# very good
