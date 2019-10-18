import numpy as np
from os import environ
from numpy.lib.stride_tricks import as_strided

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
    ret = as_strided(x, shape=(window_size, n_hops),
                          strides=(x.itemsize, hop_size * x.itemsize),
                          writeable=False)
    return ret

def get_env(name,default=None,conv=None,check_if_none=False):
    try:
        ret=environ[name]
    except KeyError:
        ret=default
        if check_if_none and ret is None:
                raise Exception("Specify " + name + ".")
        return ret
    if conv is not None:
        return conv(ret)
    return ret
