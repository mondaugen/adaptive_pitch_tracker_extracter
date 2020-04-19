import numpy as np
import common
import matplotlib.pyplot as plt

X0=common.get_env("X0",default=-10,conv=int)
Y0=common.get_env("Y0",default=10,conv=int)
X1=common.get_env("X1",default=10,conv=int)
Y1=common.get_env("Y1",default=20,conv=int)

x=np.fromfile("/tmp/dspm_2dline_s48q16_lookup_vs48q16_input.s48q16",dtype='int64')/(2**16)
y_py=(x-X0)*(Y1-Y0)/(X1-X0)+Y0
y_c=np.fromfile("/tmp/dspm_2dline_s48q16_lookup_vs48q16_output.s48q16",dtype='int64')/(2**16)

plt.plot(x,y_py,label='python')
plt.plot(x,y_c,label='c')
plt.plot(x,np.abs(y_py-y_c),label='error')
plt.legend()
plt.show()
