import numpy as np
from scipy import signal, interpolate
import matplotlib.pyplot as plt
import filters
from matplotlib.widgets import Slider

N=1000
def get_interp(a0,a1):
    return interpolate.interp1d(
        [0,N],
        np.concatenate((a0[:,None],a1[:,None]),axis=1))(np.arange(N))

def build_filter_resp(f,q,b=np.array([1]),a=np.array([1])):
    if len(f) == 0:
        return (b,a)
    b_,a_=signal.iirpeak(f[0],q[0],fs=1)
    b=np.convolve(b,b_)
    a=np.convolve(a,a_)
    return build_filter_resp(f[1:],q[1:],b,a)

f0=[0.05,0.1,0.15]
f1=[0.25,.35,.45]
q0=[10,30,100]
q1=[1,3,10]
w=np.arange(0,0.5,0.001)
b0,a0=build_filter_resp(f0,q0)
r0,c0=filters.b_a_to_r_c(b0,a0)
b1,a1=build_filter_resp(f1,q1)
r1,c1=filters.b_a_to_r_c(b1,a1)
r=get_interp(r0,r1)
c=get_interp(c0,c1)

def get_w_h(idx):
    idx=int(idx)
    r_=r[:,idx]
    c_=c[:,idx]
    if np.any(np.abs(r_)>1):
        print("unstable")
    b_,a_=filters.c_r_to_b_a(c_,r_)
    w_,h_=signal.freqz(b_,a_,worN=w,fs=1)
    return (w_,h_)

w_,h_=get_w_h(0)
plot_lines,=plt.plot(w_,20*np.log10(np.abs(h_)))
plot_lines._axes.set_ylim(-80,0)

axcolor = 'lightgoldenrodyellow'
ax_idx = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
idx_slider = Slider(ax_idx,'Index',0,N-1,valinit=0,valstep=1)

def update(val):
    idx=idx_slider.val
    w_,h_=get_w_h(idx)
    plot_lines.set_ydata(20*np.log10(np.abs(h_)))

idx_slider.on_changed(update)

plt.show()
