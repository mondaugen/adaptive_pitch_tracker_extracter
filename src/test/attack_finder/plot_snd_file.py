import numpy as np
import matplotlib.pyplot as plt
import common

path=common.get_env('FILEPATH',default="/tmp/attacks_from_sd_test_in.f32")
x=np.fromfile(path,dtype='float32')
n=np.arange(len(x))

plt.plot(n,x)
plt.show()

