import numpy as np
import matplotlib.pyplot as plt
from os import environ

vars_file='/tmp/dump_pre_and_post_transform.vars'
if 'vars_file' in environ.keys():
    vars_file=environ['vars_file']

with open(vars_file,'r') as f:
    for line in f.readlines():
        exec(line)

print("W=%d" % (W,))
print("z_input0=%d" % (z_input0,))
print("z_inputH=%d" % (z_inputH,))
print("z_outputH=%d" % (z_outputH,))

# because we use real data the length is W//2+1
Wz=W//2+1

# What the algorithm gave us
z_i0=np.fromfile("/tmp/%d.z64" % (z_input0,),dtype="complex64").reshape((-1,Wz))
z_iH=np.fromfile("/tmp/%d.z64" % (z_inputH,),dtype="complex64").reshape((-1,Wz))
print(z_iH.shape)

# What it should be (assuming np.fft is correct)
r_i0=np.fromfile("/tmp/%d.f32" % (z_input0,),dtype="float32").reshape((-1,W))
r_iH=np.fromfile("/tmp/%d.f32" % (z_inputH,),dtype="float32").reshape((-1,W))
z_i0_np=np.fft.rfft(r_i0)
z_iH_np=np.fft.rfft(r_iH)

print("max z_iH abs error", np.max(np.abs(z_iH-z_iH_np)))
print("max z_i0 abs error", np.max(np.abs(z_i0-z_i0_np)))

print("z_iH has Nan?",np.any(np.isnan(z_iH)))
print("z_i0 has Nan?",np.any(np.isnan(z_i0)))

print("z_iH has 0?",np.any(0 == z_iH))
print("z_i0 has 0?",np.any(0 == z_i0))


