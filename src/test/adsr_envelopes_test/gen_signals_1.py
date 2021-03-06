import numpy as np

N=5000
gate=np.zeros(N,dtype='float32')
gate[255:511]=1
gate[1024:1279]=1
gate[3000:4000]=1
attack_duration=np.linspace(1,1000,N).astype('uint32')
decay_duration=np.linspace(1,1000,N).astype('uint32')
#sustain_level=np.ones(N,dtype='float32')*0.5
sustain_level=np.linspace(0,1,N).astype('float32')
release_duration=np.linspace(1,1000,N).astype('uint32')

gate.tofile("/tmp/gate.f32")
attack_duration.tofile("/tmp/attack_duration.u32")
decay_duration.tofile("/tmp/decay_duration.u32")
sustain_level.tofile("/tmp/sustain_level.f32")
release_duration.tofile("/tmp/release_duration.u32")
