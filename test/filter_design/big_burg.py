import filters
import numpy as np
from scipy import signal

b=[1]
P=300
R=np.random.uniform(-1,1,P)
a=filters.r_to_A(R)[-1,:]
N=6000
x=np.random.standard_normal(N)
x=signal.lfilter(b,a,x)

# estimate gamma
#gamma,eplus,eminus=filters.burg2(x,P)
gamma,eplus,eminus=filters.burg(x,P)
a_=filters.r_to_A(gamma)[-1,:]
x_=eplus
x_=signal.lfilter([1],a_,x_)
error_a=signal.lfilter(a_,[1],x)

#print('R',R)
#print('gamma',gamma)
print(np.mean(np.abs(R-gamma)))
