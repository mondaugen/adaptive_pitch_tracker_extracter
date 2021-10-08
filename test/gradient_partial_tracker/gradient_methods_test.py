import numpy as np
import matplotlib.pyplot as plt
from some_ft import dirichlet, dirichlet_dk
from some_sig import rectangular, sinusoid_about_0
from dftdk import dft_dk, dft_bin, ps_dk, gradient_ascent_step, newton_ascent_step, self_adjusting_ga_step

j=complex('j')
N=4096
n=np.arange(N)
W=2047
n=np.arange(N)
# The true peak position
k_star=1234
x=rectangular(n,W,N)*sinusoid_about_0(k_star,N)+np.random.standard_normal(N)*0.003
X=dft_dk(N,p=0)(x)
plt.plot(n,np.real(X),label="Fourier transform of x")
mu_ga=1
mu_newton=0.1 # ?
mu_saga=0.3 # ?
n_steps=10

# now try and find the peak
def do_gradient_method(k0,n_steps,grad_meth):
    k=np.zeros(n_steps+1)
    k[0]=k0
    for step in range(n_steps):
        k0 = grad_meth(x,k0)
        k[step+1]=k0
    return k

k0=1236
k_ga=do_gradient_method(k0,n_steps,lambda x,k0: gradient_ascent_step(x,k0,mu_ga))
k_newton=do_gradient_method(k0,n_steps,lambda x,k0: newton_ascent_step(x,k0))
plt.plot(k_ga,np.real(dft_bin(x,k_ga)),label='gradient ascent')
plt.plot(k_newton,np.real(dft_bin(x,k_newton)),label='newton')

k_saga=np.zeros(n_steps+1)
saga_steps=np.zeros(n_steps+1)
k_saga[0]=k0
step_saga=1
saga_steps[0] = step_saga
for i in range(n_steps):
    step_saga,k0 = self_adjusting_ga_step(x,k0,step_saga,1,mu_saga)
    k_saga[i+1] = k0
    saga_steps[i+1] = step_saga
print(saga_steps)
plt.plot(k_saga,np.real(dft_bin(x,k_saga)),label='saga')
    
plt.legend()

plt.xlim(1230,1240)
plt.show()
