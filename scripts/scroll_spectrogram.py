"""

Load in sounfile
For each hop, plot the spectrum. This constitutes a frame.
Concatenate the frames to make a movie.

"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import signal

sr=16000
x=np.fromfile('/tmp/snd.f64')
H=256
N_FFT=4096
N_W=1024
N_x=len(x)
w=signal.get_window('hann',N_W)
W=np.sum(w)
x_w=np.zeros(N_FFT)

def data_gen(t=0):
    cnt = 0
    while cnt < (len(x)-N_W):
        x_w[(N_FFT-N_W)//2:(N_FFT-N_W)//2+N_W]=x[cnt:cnt+N_W]*w
        X_=20*np.log10(np.abs(np.fft.fft(x_w)/W))
        t = (cnt + N_W) / sr
        cnt += H
        yield (t, X_)

fig, ax = plt.subplots()
ax.grid()
xdata=np.arange(N_FFT/2)*sr/N_FFT
ydata=[]

def init():
    pass

def run(data):
    # update the data
    t, y = data
    ax.clear()
    ax.set_ylim(-160, 3)
    ax.set_xlim(0, sr/2)
    ax.plot(xdata,y[:N_FFT//2],label="%f" % (t,))
    ax.legend()

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15)

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=200,
                              repeat=False, init_func=init)

#ani.save('/tmp/spec.mp4',writer=writer)
plt.show()
