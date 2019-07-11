import numpy as np
from os import environ

def normalize(x):
    x-=np.mean(x)
    x/=np.max(np.abs(x))
    return x

def frame(x,hop_size,window_size):
    """
    A "frame" function independent of librosa.
    x must be one dimensional.
    windows are output in the colums, so that the whole resulting matrix can be
    left-multiplied to perform a transform.
    """
    n_hops=(len(x)-window_size)//hop_size
    n=np.add.outer(np.arange(window_size),np.arange(n_hops)*hop_size)
    ret=x[n]
    return ret

def get_env(name,default,conv=None):
    try:
        ret=environ[name]
    except KeyError:
        return default
    if conv is not None:
        return conv(ret)
