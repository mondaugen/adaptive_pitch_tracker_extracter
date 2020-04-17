# Compute what fixed point precision we need to acheive a certain speed to
# within certain accuracy

import numpy as np
import matplotlib.pyplot as plt

n=np.arange(1,32)
rate=(2**(1/1200)-1)*(2**n)
rounded_rate=np.round(rate)
err = np.log10(np.abs(rounded_rate/rate-1))
plt.plot(n,err)
plt.show()
