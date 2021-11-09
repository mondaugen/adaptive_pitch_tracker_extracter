import matplotlib.pyplot as plt
from _common import check_dv_dft
from some_sig import multi_mod_sum_of_cos
from some_ft import multi_mod_sum_of_cos_dft_k, multi_mod_sum_of_cos_dft_dk, multi_mod_sum_of_cos_dft_k
from functools import partial
import numpy as np

p=(1+np.arange(10))
B=1./p
V=0.05*p
check_dv_dft(
    N=4096,
    W=2047,
    L=2048,
    v=V,
    x_td=lambda v,A,L,W,N: multi_mod_sum_of_cos(B,v,A,L,W,N),
    x_dft=lambda k0,k,A,L,W,N: multi_mod_sum_of_cos_dft_k(B,k0,k,A,L,W,N),
    x_ddft=lambda k0,k,A,L,W,N: multi_mod_sum_of_cos_dft_dk(B,k0,k,A,L,W,N)
)
plt.show()
