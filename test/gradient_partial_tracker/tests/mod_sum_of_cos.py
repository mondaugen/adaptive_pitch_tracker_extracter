from some_sig import mod_sum_of_cos
from some_ft import mod_sum_of_cos_dft, normalize_sum_of_cos_A
import numpy as np

def compare(
    v=0.25,
    N=8,
    L=8,
    W=7,
    A=[0.5,0.5],
    normalize=True
):
    if normalize:
        A=normalize_sum_of_cos_A(A,L,W,N)
    k=np.arange(N)
    w=mod_sum_of_cos(v,A,L,W,N)
    w_ft=np.fft.fft(w)
    w_ft_th=mod_sum_of_cos_dft(v,k,A,L,W,N)

    print("w_ft")
    print(np.real(w_ft))
    print(np.imag(w_ft))
    print("w_ft_th")
    print(np.real(w_ft_th))
    print(np.imag(w_ft_th))

compare()
compare(v=0.33)
compare(N=16)
