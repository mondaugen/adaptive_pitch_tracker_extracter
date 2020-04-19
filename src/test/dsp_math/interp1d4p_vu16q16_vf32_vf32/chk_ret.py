import numpy as np
import matplotlib.pyplot as plt
import vectorized_cubic

x=np.fromfile("/tmp/x.f32",dtype="float32")
y=np.fromfile("/tmp/y.f32",dtype="float32")
xi=np.fromfile("/tmp/xi.u16q16",dtype="uint32")
xr=np.fromfile("/tmp/xr.u32",dtype="uint32")

# python implementation
yi=vectorized_cubic.lookup4(xi/(2**16),y)
err=yi-y[xr]
err=20*np.log10(np.abs(err))

# c implementation
yi_c=np.fromfile("/tmp/yi_c.f32",dtype="float32")
err_c=yi_c-y[xr]
err_c=20*np.log10(np.abs(err_c))
err_c[err_c<-100]=-100

fig,ax=plt.subplots(2,1)
ax[0].plot(x*(2**16),y,label='original')
ax[0].plot(xi,yi,label='interpolated')
ax[0].plot(xi,yi_c,label='interpolated c')
ax[0].legend()
ax[1].plot(xr,err,label='err dB')
ax[1].plot(xr,err_c,label='err c dB')
ax[1].legend()

plt.show()

