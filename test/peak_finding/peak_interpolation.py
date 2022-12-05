import cubic_sinusoid_synth as cbs
import numpy as np
import common

# phases
Fs=16e3
T=1
f0=100
fT=100
P=50
H=256
n_bp=int(np.ceil(Fs*T/H))
omega=np.multiply.outer(
    2*np.pi/Fs*np.linspace(f0,fT,n_bp),
    np.arange(P)+1
)
theta=np.zeros(P)
ph=cbs.quadratic_phase_poly_interp(H,theta,omega)

# amplitudes
a0=0
aT=-60
fall_off=-12
p_amps=fall_off*np.log2(np.arange(P)+1)
print(p_amps)
a=np.add.outer(
    np.linspace(a0,aT,n_bp),
    p_amps
)
A=cbs.linear_amplitude_interp(H,a)
x=np.sum(np.cos(ph)*np.power(10.,A/20.),axis=0)
x=common.normalize(x)
x.tofile('/tmp/out.f64')

