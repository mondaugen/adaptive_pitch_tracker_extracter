import numpy as np
from os import environ
from numpy.lib.stride_tricks import as_strided
import uuid
import matplotlib.pyplot as plt

def normalize(x,subtract_mean=True):
    if subtract_mean:
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

def frame_times(len_x,H,W):
    """ Get the sample times of the beginning of each hop when framing """
    # Right now a hack is used some day maybe we'll actually do a calculation
    fake_x=np.ones(len_x)
    fake_x_framed=frame(fake_x,H,W)
    return np.arange(fake_x_framed.shape[1])*H

def ola_shorten(x,H):
    """
    Make a signal of length H by taking consecutive sections of x of length H
    and summing them into an output of length H. If there exists a section that
    would exceed the end of x, the values outside of x are taken to be 0.
    """
    W=len(x)
    x_ext=np.concatenate((x,np.zeros(H-W%H)))
    r=np.zeros(H)
    for h in np.arange(0,len(x),H):
        r+=x_ext[h:h+H]
    return r

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

def next_pow_2(n):
    x = 1
    while x < n:
        x = x << 1
    return x

def mktemp(template='/tmp/%s%s',suf=''):
    """ Make a temporary file. Not secure. """
    return template % (str(uuid.uuid4()),suf)

def logic_plot(x,y,ax=None,**kwargs):
    len_x_=len(x)*2-1
    x_=np.zeros(len_x_,dtype=x.dtype)
    y_=np.zeros(len_x_,dtype=y.dtype)
    x_[::2]=x
    x_[1::2]=x[1:]
    y_[::2]=y
    y_[1::2]=y[:-1]
    if ax is None:
        plt.plot(x_,y_,**kwargs)
    else:
        ax.plot(x_,y_,**kwargs)

def region_plot(start,end,height=1,level=0,ax=None,**kwargs):
    x=[start,start+0.5,end-0.5,end]
    y=[level,level+height,level+height,level]
    if ax is None:
        plt.plot(x,y,**kwargs)
    else:
        ax.plot(x,y,**kwargs)

# plot an arch showing where the analysis window is placed
def plot_arch(ax,s,l,h,b=0,color='black',print_h=None,t=0):
    # t is taper
    ax.plot([s,s+t],[b,b+h],color=color)
    ax.plot([s+t,s+l-t],[b+h,b+h],color=color)
    ax.plot([s+l-t,s+l],[b+h,b],color=color)
    if print_h is not None:
        ax.text(s,b+h, print_h % (h,))
