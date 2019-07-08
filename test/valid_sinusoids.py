# Plot a random sinusoid that is a feasible solution

import numpy as np
import matplotlib.pyplot as plt

show_plot=False

N=10000

j=complex('j')
# frequency values
w0=2*np.pi*440./16000.*np.power(2,np.random.uniform(-12,12)/12)
wmin=w0*np.power(2,-7/12)
wmax=w0*np.power(2,7/12)
dwmin=-3/N*10
dwmax=3/N*10
# amplitude values
amin=np.log(1e-3)
amax=np.log(1)
damin=np.log(1e-8)/N*30
damax=np.log(1e8)/N*30
# order of polynomial
p_ord=5

n=np.arange(N)

dw_vals=w0*(np.power(2,np.random.uniform(dwmin,dwmax,(N,)))-1)
w_vals=np.random.uniform(wmin,wmax)+np.cumsum(dw_vals)
phi_vals=np.cumsum(w_vals)+np.random.uniform(-np.pi,np.pi)
print(phi_vals)
da_vals=np.random.uniform(damin,damax,size=(N,))
a_vals=np.cumsum(da_vals)+np.random.uniform(amin,amax)
print(a_vals)

p_a=np.polyfit(n,a_vals,p_ord)
p_phi=np.polyfit(n,phi_vals,p_ord)

y=np.exp(a_vals+j*phi_vals)
y_=np.exp(np.polyval(p_a,n)+j*np.polyval(p_phi,n))
y_mod=np.exp(np.polyval(p_a,n))

plt.plot(n,np.real(y_))
plt.plot(n,np.real(y))
plt.plot(n,y_mod)

out=y#np.concatenate((y,np.real(y_)))
out-=np.mean(out)
out/=np.max(np.abs(out))
with open('/tmp/exp.f64','a') as fd:
    out.tofile(fd)

if show_plot:
    plt.show()
