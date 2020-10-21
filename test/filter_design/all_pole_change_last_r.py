# See what happens to the filter poles when you change the last reflection coefficient

import filters
import numpy as np
import matplotlib.pyplot as plt

def adjust_r(A,last_r_rad,replace=False):
    # find reflection coefficients
    R=filters.a_to_r(A)
    R[-1]/=np.abs(R[-1])
    if replace:
        R[-1]=last_r_rad
    else:
        R[-1]*=last_r_rad
    Az_new=filters.r_to_A(R)[-1,:]
    return Az_new

def plot_poles_from_poly(A,label):
    poles_new=1/np.roots(A[::-1])
    poles_real_new=np.real(poles_new)
    poles_imag_new=np.imag(poles_new)
    plt.plot(poles_real_new,poles_imag_new,'.',label=label)

# radius and angle (normalized) of poles
poles_rad_angle=[(1-1e-5,v) for v in np.arange(0.01,0.16,0.01)]

poles_complex=[r*np.exp(2*np.pi*v*1j) for r,v in poles_rad_angle]

# build up filter
Az=1
for pole in poles_complex:
    Az=np.convolve(Az,[1,-pole])

# R[-1] forced radius
R_1_rads = np.interp(np.linspace(0,1,100),[0,1],[1-1e-2,1-1e-6])
R_1_replace=False


# show poles
plot_poles_from_poly(Az,'original')

# plot unit circle

theta = np.linspace(-np.pi, np.pi, 201)
plt.plot(np.sin(theta), np.cos(theta), color = 'gray', linewidth=0.2)

#for n,new_r in enumerate(R_1_rads):
#    Az_new=adjust_r(Az,new_r,replace=R_1_replace)
#    plot_poles_from_poly(Az_new,label='%d'%(n,))

plt.grid()
plt.legend()
plt.show()
