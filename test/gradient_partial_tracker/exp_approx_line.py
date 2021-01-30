import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

# Approximate a line with an exponential

# The range of the resulting line
N=1024
# maximum error tolerated
err_max=1e-3/N
# We want to find the point at which the error function exp(x) - (x+1) ==
# err_max, so we find the root of this function:
def f(x):
    return np.exp(x) - (x + 1) - err_max
# Its derivative
def df(x):
    return np.exp(x) - 1
sol0 = optimize.root_scalar(f, method='newton', fprime=df, x0=-1)
print(sol0.root)
sol1 = optimize.root_scalar(f, method='newton', fprime=df, x0=1)
print(sol1.root)
print('exp(x0)',np.exp(sol0.root))
print('exp(x1)',np.exp(sol1.root))
print('decay rate exp((x1-x0)/N)',np.exp(-(sol1.root-sol0.root)/N))
points=np.linspace(sol0.root,sol1.root,N)
plt.plot(points,points+1,label='line')
plt.plot(points,np.exp(points),label='exp')
plt.legend()
plt.show()
