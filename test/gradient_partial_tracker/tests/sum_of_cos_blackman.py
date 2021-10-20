from some_sig import sum_of_cos
import matplotlib.pyplot as plt
import numpy as np

N=64
w=sum_of_cos(np.array([0.42,0.5,0.08]),N)
n=np.arange(N)
print(w[N//2])

plt.plot(n,w)
plt.show()

