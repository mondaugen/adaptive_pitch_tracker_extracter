# Test signals

import cubic_sinusoid_synth as cbs
import numpy as np
import common

def sum_of_sinusoids(
    # time
    Fs=16e3,
    T=1,
    H=256,
    # n partials
    P=50,
    # phase
    f0=100,
    fT=100,
    # amplitudes
    a0=0,
    aT=-60,
    fall_off=-12,
    # noise
    SNR=-60
):

    n_bp=int(np.ceil(Fs*T/H))
    omega=np.multiply.outer(
        2*np.pi/Fs*np.linspace(f0,fT,n_bp),
        np.arange(P)+1
    )
    theta=np.zeros(P)
    ph=cbs.quadratic_phase_poly_interp(H,theta,omega)
    N=ph.shape[1]

    p_amps=fall_off*np.log2(np.arange(P)+1)
    a=np.add.outer(
        np.linspace(a0,aT,n_bp),
        p_amps
    )
    A=cbs.linear_amplitude_interp(H,a)

    # sinusoids
    s=np.sum(np.cos(ph),axis=0)
    s_pow=np.sum(s**2)/N
    no_g=s_pow*(10**(SNR/20))
    no=np.random.normal(size=(N,))*(no_g**0.5)
    x=np.sum((no+np.cos(ph))*np.power(10.,A/20.),axis=0)
    x=common.normalize(x)

    return x

