# Processing of single channel numpy arrays for synthesis
import numpy as np
from scipy import signal

def cutoff_to_normalized(cutoff,sample_rate):
    """ Convert the frequency to the frequency that the filter design functions
    expect, which is 0 for DC, 1 for Nyquist """
    return cutoff/(sample_rate*0.5)

def lowpass_1pole_filter(x,cutoff,sample_rate):
    freq=cutoff_to_normalized(cutoff,sample_rate)
    b,a=signal.butter(1,freq)
    return signal.lfilter(b,a,x)

def gain_db(x,gain,sample_rate):
    gain_lin=10**(gain/20)
    return x*gain_lin

def truncate(x,length,sample_rate,fadetime=0.25,min_db=-120):
    """
    Fade is exponential, so values below min_db are changed to 0. min_db is
    also the lowest point of the fade
    """
    if len(x) < 1:
        return x
    fade_start_samples=int(np.round((length-fadetime)*sample_rate))
    if fade_start_samples < 0:
        fade_start_samples = 0
    fade_end_samples=int(np.round(length*sample_rate))
    if fade_end_samples > len(x):
        fade_end_samples = len(x)
    # Not sure if the linear interpolation function is happy if we have two
    # equal x values, so we offset the start and end of the fade by a fraction
    # of a sample
    fade_curve=np.interp(
    np.arange(len(x)),
    [0,fade_start_samples+1e-3,fade_end_samples-1e-3,len(x)],
    [0,0,min_db,min_db])
    fade_curve_lin=np.power(10,fade_curve/20)
    fade_curve_lin[fade_curve_lin < np.power(10,min_db/20)]=0
    return x*fade_curve_lin

proc_funs={
'lowpass_1pole': lowpass_1pole_filter,
'gain_db': gain_db,
'truncate': truncate
}

def parse_proc_fun(x,s):
    funname=s.split()[0]
    args=",".join(s.split()[1:])
    dict_args="dict(%s)" % (args,)
    return proc_funs[funname](x,**eval(dict_args))

class const_sr_fun_dispatcher:
    """ So you don't always have to specify the sample rate as an argument, this
    will append it for you when you parse a processing function. """
    def __init__(self,sample_rate):
        self.sample_rate=sample_rate
    def parse_proc_fun(self,x,s):
        s+=" sample_rate=%f" % (self.sample_rate,)
        return parse_proc_fun(x,s)
        
