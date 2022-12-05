import tsigs
from os import environ as env
x=tsigs.sum_of_sinusoids(
    # time
    Fs=float(env.get('Fs','16e3')),
    T=float(env.get('T','1')),
    H=int(env.get('H','256')),
    # n partials
    P=int(env.get('P','50')),
    # phase
    f0=float(env.get('f0','100')),
    fT=float(env.get('fT','100')),
    # amplitudes
    a0=float(env.get('a0','0')),
    aT=float(env.get('aT','-60')),
    fall_off=float(env.get('fall_off','-12')),
    # noise
    SNR=float(env.get('SNR','-60'))
)
OUTFILE=env.get('OUTFILE','/tmp/out.f64')
x.tofile(OUTFILE)

