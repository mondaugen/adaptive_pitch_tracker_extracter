import filters
import numpy as np
from scipy import signal

b=[1]
a=[1,-0.12,-0.456,0.6]
N=6000
P=3
x=np.zeros(N)
x[0]=1
x=signal.lfilter(b,a,x)

# estimate gamma
#gamma,eplus,eminus=filters.burg2(x,P)
gamma,eplus,eminus,_=filters.burgN(x,P)
a_=filters.r_to_A(gamma)[-1,:]
x_=eplus
x_=signal.lfilter([1],a_,x_)
error_a=signal.lfilter(a_,[1],x)

print('gamma',gamma)
print('a_',a_)
#print('absolute resynthesis error', np.mean(np.abs(x-x_)))
print('eplus',eplus[:10])
print('error_a',error_a[:10])
print('eminus',eminus[:20])
print('x',x[:20])
print('x_',x_[:20])
