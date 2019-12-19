import numpy as np
import lfo

SR=16000
f0=0.5/SR
N=5000000
gate = lfo.chirp(N,f0,f0,phase=-0.5).squarewave().astype('float32')
attack_duration=np.linspace(100,10000,N,dtype='uint32')
decay_duration=np.linspace(10000,100,N,dtype='uint32')
sustain_level=lfo.chirp(N,.123*f0,f0,x_min=0.25,phase=-0.25).sinusoid().astype('float32')
#sustain_level=np.linspace(0,1,N).astype('float32')
release_duration=lfo.chirp(N,f0,.321*f0,x_min=100,x_max=100000,phase=-0.25).sinusoid().astype('uint32')

gate.tofile("/tmp/gate.f32")
attack_duration.tofile("/tmp/attack_duration.u32")
decay_duration.tofile("/tmp/decay_duration.u32")
sustain_level.tofile("/tmp/sustain_level.f32")
release_duration.tofile("/tmp/release_duration.u32")

