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


N=2048
W=1023
L=W+1
k0_sig=N*0.1
k0_anl=N*np.power(0.1,1/12) # can specify in semitones
# partial multipliers
p=(1+np.arange(10))
# amplitude scalars
B=1./p
A=normalize_sum_of_cos_A([0.42,0.5,0.08],L,W,N)
dX_dk=multi_mod_sum_of_cos_dft_dk(B,k0_sig*p,k0_anl*p,A,L,W,N)
avg_dX_dk=dX_dk.mean()

print(avg_dX_dk)

# TODO: Find the same value simply by finding the inner product with a signal
# that is the sum of all the harmonics
