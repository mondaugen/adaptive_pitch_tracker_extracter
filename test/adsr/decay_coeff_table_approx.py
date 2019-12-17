# how well does linear interpolation of a tabulated pow(a,1/n) work?

import numpy as np
import common
import matplotlib.pyplot as plt

decay_min_dB=-60
decay_min_A=np.power(10,decay_min_dB/20)
decay_coeff_tab_max=5*48000
N=common.next_pow_2(decay_coeff_tab_max)
x=np.arange(N)
y=np.power(decay_min_A,1./x)
len_tab=int(np.log2(N))
x_=np.power(2,np.arange(len_tab+1))
y_=np.power(decay_min_A,1./x_)

plt.plot(x,y,label='true')
plt.plot(x_,y_,label='linear interpolation')
plt.legend()

plt.show()

