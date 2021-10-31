from some_ft import sum_of_cos_dft
from some_sig import sum_of_cos
from scipy import signal
import numpy as np
#w = sum_of_cos([1],8,7,8)
L=8
W=7
N=8
A=[0.5,0.5]
w = sum_of_cos(A,L,W,N)
print(w)
print("something interesting:")
w2 = np.concatenate((w[-(N//2):],w[:N//2]))
print("w[-(N//2):] ++ w[:N//2]:",w2)
print("signal.get_window('hann',8)")
print(signal.get_window("hann",8))
print("So it depends on what you define as the 'length of the window!'")
# The window is made of a sum of cosines windowed by a boxcar of length {N} but
# the period of those cosines is {L}. scipy.signal uses this period value as the
# window length and considers the 0 as part of the window, but if you are
# interested in zero-padding the resulting signal, it makes more sense to
# consider only the length of contiguous non-zero values (inside the boxcar) as
# the length (as we have done here with {W}).

fw=np.fft.fft(w)
k=np.arange(N)
fw_th=sum_of_cos_dft(k,A,L,W,N)
print("fw:",np.real(fw))
print("fw_th:",fw_th)

w_zp=np.concatenate((w[:N//2],np.zeros(N),w[-(N//2):]))
fw_zp=np.fft.fft(w_zp)
k=np.arange(2*N)
fw_zp_th=sum_of_cos_dft(k,A,L,W,2*N)
print("fw:",np.real(fw_zp))
print("fw_th:",fw_zp_th)
# looking good
