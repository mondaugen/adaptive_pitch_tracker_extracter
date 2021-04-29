import numpy as np
from polyeval import polyeval

P=4
K=3
N=10
pcoefs=np.random.standard_normal((P,K))
z=np.linspace(-1,1,N)
r_polyeval=polyeval(pcoefs,z)
r_numpy=np.zeros((K,N))
for k,pcoefs_ in enumerate(pcoefs.T):
    r_numpy[k,:]=np.poly1d(pcoefs_)(z)

print(np.abs(r_polyeval - r_numpy).max())

print("r_numpy") 
print(r_numpy) 
print("r_polyeval") 
print(r_polyeval) 
