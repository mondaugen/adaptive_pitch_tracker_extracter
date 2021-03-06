# Tools for windowing signals, or getting parts of signals.
import numpy as np
from scipy import signal
from math import ceil

def last_slice_index(s):
    ''' Get the last value a slice object can give (why is this not builtin?) '''
    ret = (ceil((s.stop-s.start)/s.step)-1)*s.step+s.start
    if ret < s.start:
        # This occurs if s.stop is <= s.start
        return None
    return ret

class signal_windower:
    def __init__(self,x,oob_noise_power=-100):
        """
        x is the array you want to index
        oob_noise_power is the power in dB of the out of bounds noise
        """
        if len(x) < 1:
            raise ValueError('x must have positive length')
        self.x=x
        # The random key for generating out-of-bound values on the left-hand
        # side (so they are always the same)
        self._left_random_key=np.random.randint(1e6)
        # The random key for generating out-of-bound values on the right-hand
        # side (so they are always the same)
        self._right_random_key=np.random.randint(1e6)
        # This is because to get a desired power of normally distributed values,
        # you multiply by the square root of the desired power.
        self.oob_noise_scalar=10**(oob_noise_power/40)
    def __getitem__(self,key):
        if type(key) == type(int()):
            if key < 0:
                lh_rs=np.random.RandomState(self._left_random_key)
                lh_vals=lh_rs.standard_normal(-key)
                return lh_vals[-key-1]*self.oob_noise_scalar
            if key >= len(self.x):
                rh_rs=np.random.RandomState(self._right_random_key)
                rh_vals=rh_rs.standard_normal(key-len(self.x)+1)
                return rh_vals[key-len(self.x)]*self.oob_noise_scalar
            return self.x[key]
        if type(key) == type(slice(0,0,1)):
            min_key=key.start
            max_key=last_slice_index(key)
            ret=self.x
            slice_start=key.start
            slice_stop=key.stop
            if min_key < 0:
                lh_rs=np.random.RandomState(self._left_random_key)
                lh_vals=lh_rs.standard_normal(-min_key)[::-1]*self.oob_noise_scalar
                ret=np.concatenate((lh_vals,ret))
                slice_start=0
                slice_stop=slice_stop-min_key
            if max_key >= len(self.x):
                rh_rs=np.random.RandomState(self._right_random_key)
                rh_vals=rh_rs.standard_normal(max_key-len(self.x)+1)*self.oob_noise_scalar
                ret=np.concatenate((ret,rh_vals))
            return ret[slice(slice_start,slice_stop,key.step)]

class taper_window_applier:
    """ Makes it easy to apply a window to the beginning or end of a signal """
    def __init__(self,window_type='hann',table_length=4096):
        self.window=signal.get_window(window_type,table_length)
        self.table_length=table_length
    def get_taper(self,N,where='beginning'):
        """
        Get a taper function of length N, which will be values going from 0 to
        the middle value of the window.
        if where is 'beginning', we fade in, if 'end' we fade out.
        if where is unrecognized, we use 'beginning'
        """
        r=np.interp(
            np.linspace(0,self.table_length//2,N),
            np.arange(self.table_length//2+1),
            self.window[:self.table_length//2+1])
        if where == 'end':
            return r[::-1]
        return r

def windowed_lookup_default_fill(W):
    return np.random.standard_normal(W)*1e-6

class windowed_lookup:
    """
    A port of windowed_lookup_f32.c
    """
    def __init__(self,
        # signal to look up from
        x,
        # maximum size of window used for lookup
        W,
        # function used to fill out of bounds values
        fill_func=windowed_lookup_default_fill):
        self.W=W
        self.x=x
        self.len_x=len(x)
        self.signal_start=np.concatenate((fill_func(W).astype(x.dtype),x[:W]))
        self.signal_end=np.concatenate((x[-W:],fill_func(W).astype(x.dtype)))
    def access(self,index):
        if index < -self.W:
            index = self.W
        if index > self.len_x:
            index = self.len_x
        if index < 0:
            ret = self.signal_start[self.W+index:2*self.W+index]
        elif index > (self.len_x - self.W):
            ret = self.signal_end[self.W-(self.len_x-index):2*self.W-(self.len_x-index)]
        else:
            ret = self.x[index:index+self.W]
        return ret
