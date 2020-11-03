import numpy as np
from scipy import signal
from spectral_difference import local_max, local_min

def hold_values(x,d=0,v=None):
    # run the non-linear filter:
    # y[n] = y[n-1] * (v[n] == d) + x[n] * (v[n] != d)
    # which holds the last value over values that are d, or otherwise introduces
    # a new value.
    # if v is none, it defaults to x
    # TODO: write this routine in C
    if v is None:
        v = x
    y=np.zeros_like(x)
    yn_1=0
    for n in range(len(x)):
        y[n] = yn_1 * (v[n] == d) + x[n] * (v[n] != d)
        yn_1 = y[n]
    return y

def find_peaks(x,K=1,T=float('-inf')):
    # find peaks that are K times greater than the nearest minimum
    # T is the absolute threshold, underwhich peaks are not output
    max_i=local_max(x,one_sided_max='none')
    min_i=local_min(x,one_sided_max='none')
    min_maxs=np.zeros_like(x)
    min_maxs[max_i]=1
    # disallow maxima on first or last points
    min_maxs[0] = 0
    min_maxs[-1] = 0
    min_maxs[min_i]=-1
    min_maxs_vals=np.zeros_like(x)
    min_maxs_vals[max_i]=x[max_i]
    min_maxs_vals[min_i]=x[min_i]
    ext_fhold=hold_values(min_maxs)
    ext_bhold=hold_values(min_maxs[::-1])[::-1]
    ext_val_fhold=hold_values(min_maxs_vals,d=0,v=min_maxs)
    ext_val_bhold=hold_values(min_maxs_vals[::-1],d=0,v=min_maxs[::-1])[::-1]
    max_over_min=np.zeros_like(x)
    # get boolean where both maxima, otherwise minima
    mask = (ext_fhold == ext_bhold) & (ext_fhold == 1)
    max_over_min[mask] = ext_val_fhold[mask]
    # get extended minima
    ext_min=(((ext_fhold == -1) * ext_val_fhold + (ext_bhold == -1) * ext_val_bhold)
        * np.power(0.5,ext_fhold == ext_bhold)) # divide by 2 if both were added
    mask = ~mask
    max_over_min[mask] = ext_min[mask]
    # from this, maxima can be filtered by comparing with adjacent values to see
    # if they are K times greater
    peaks=local_max(max_over_min,one_sided_max='none',K=K)
    peaks=peaks[max_over_min[peaks]>T]
    return peaks,max_over_min,max_i,min_i,ext_val_fhold,ext_val_bhold
