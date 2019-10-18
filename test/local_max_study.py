import numpy as np
import matplotlib.pyplot as plt

def local_max(x):
    """ for matrices, finds the local maxima within the columns """
    gtr=np.concatenate((x[:-1,:]>=x[1:,:],np.zeros((1,x.shape[1]),dtype='bool')),axis=0)
    gtl=np.concatenate((np.zeros((1,x.shape[1]),dtype='bool'),x[1:,:]>x[:-1,:]),axis=0)
    return np.where(gtr&gtl)

R=10
C=3

x=np.arange(R)
y=np.random.standard_normal((R,C))
r_max,c_max=local_max(y)
c_max_sort_n=np.argsort(c_max)
r_max=r_max[c_max_sort_n]
c_max=c_max[c_max_sort_n]

for n in range(C):
    plt.plot(x,y[:,n])

plt.plot(r_max,y[r_max,c_max],'.')

plt.show()

