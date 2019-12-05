import numpy as np
import matplotlib.pyplot as plt
import common

IN_FILE=common.get_env('IN_FILE',default='/tmp/out.f64')
SR=common.get_env('SR',default=16000,conv=int)

x=np.fromfile(IN_FILE,dtype='float64')
N=len(x)
n=np.arange(N)
t=n/SR

plt.plot(t,x)
plt.show()
