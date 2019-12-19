import numpy as np
import matplotlib.pyplot as plt

x=np.fromfile("/tmp/adsr_envelope.f32",dtype='float32')
n=np.arange(len(x))
print(np.max(x))
#nan_idx=np.where(np.isnan(x))[0][0]
#nan_idx=300
#print(nan_idx)
#x=x[:nan_idx-1]
#n=n[:nan_idx-1]
plt.plot(n,x)
plt.show()
