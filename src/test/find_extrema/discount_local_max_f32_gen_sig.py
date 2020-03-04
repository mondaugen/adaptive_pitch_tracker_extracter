import numpy as np
import matplotlib.pyplot as plt
import lfo
import common

DO_PLOT=common.get_env('DO_PLOT',default=False,conv=bool)

N=5000
n=np.arange(N)
x=lfo.chirp(N,0.01,0.02).sawtooth() 
x=x*x
x*=lfo.chirp(N,0.001,0.001).sinusoid()
x*=np.random.standard_normal(N)
x=common.normalize(x)*.99
x=np.abs(x)
x.astype('float32').tofile('/tmp/noi.f32')

if DO_PLOT:
    plt.plot(n,x)
    plt.show()
