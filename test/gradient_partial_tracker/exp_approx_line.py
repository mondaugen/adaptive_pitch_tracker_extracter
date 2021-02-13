import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

# Approximate a line with an exponential

def exp_approx_line_coeffs(a,b,err_max):
    """
    Find the coefficients of a function exp(alpha*x + beta) + gamma that
    approximates a line function starting at 0, and has value y=b at x=a with
    the maximum error err_max
    """
    # We want to find the point at which the error function exp(x) - (x+1) ==
    # err_max, so we find the root of this function:
    def f(x):
        return np.exp(x) - (x + 1) - err_max
    # Its derivative
    def df(x):
        return np.exp(x) - 1
    sol0 = optimize.root_scalar(f, method='newton', fprime=df, x0=-1)
    sol1 = optimize.root_scalar(f, method='newton', fprime=df, x0=1)
    s0 = sol0.root
    s1 = sol1.root
    alpha=(s1-s0)/a
    beta=s0+np.log(b/(np.exp(s1)-np.exp(s0)))
    gamma=-np.exp(s0)/(np.exp(s1)-np.exp(s0))
    return (alpha,beta,gamma,s0,s1)

# The range of the resulting line
N=1024
# maximum error tolerated
err_max=1e-4/N
alpha,beta,gamma,s0,s1=exp_approx_line_coeffs(N,1,err_max)
print(s0)
print(s1)
print('exp(x0)',np.exp(s0))
print('exp(x1)',np.exp(s1))
print('decay rate exp((x1-x0)/N)',np.exp(-(s1-s0)/N))
points=np.linspace(s0,s1,N)
fig,ax=plt.subplots(2,1)
ax[0].plot(points,points+1,label='line')
ax[0].plot(points,np.exp(points),label='exp')
ax[0].legend()
ax[0].set_title('line vs exp')
# now show how it performs in the sitation where you are using it to approximate
# the line described by a and b
n=np.arange(N)
y=n/N
ax[1].plot(n,y,label='line')
ax[1].plot(n,np.exp(alpha*n+beta)+gamma,label='exp')
ax[1].legend()
ax[1].set_title('offset line vs offset exponential')
print('max error for offset functions:',np.max(np.abs(y-(np.exp(alpha*n+beta)+gamma))))
print('number of bits for decay coefficient:',np.log2(1-np.exp(-(s1-s0)/N)))
plt.tight_layout()
plt.show()
