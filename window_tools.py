# Tools for windowing signals, or getting parts of signals.
import numpy as np
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

                
            
