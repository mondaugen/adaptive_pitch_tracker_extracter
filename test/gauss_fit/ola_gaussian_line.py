# See if overlapping and adding gaussians makes a constant 1

import numpy as np
import matplotlib.pyplot as plt

def g(x,K=10,a=1):
    k=np.arange(-K,K+1)*a
    G=np.exp(-1*np.power(np.add.outer(-1*k*np.sqrt(-1*np.log(0.5)),x),2))
    return G.sum(axis=0)
    
x=np.arange(-100,100,0.1)
print('K=10',g(0))
print('K=100',g(0,K=100))
print('K=1000',g(0,K=1000))
plt.plot(x,g(x,K=100,a=1))
plt.show()
