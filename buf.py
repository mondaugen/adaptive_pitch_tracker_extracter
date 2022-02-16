# Datatypes storing sampled signals

import numpy as np

class inf_buf:
    """
    An infinitely long buffer. Reads outside of some initialized region return a
    fill value (default 0).
    """
    def __init__(self,s,start=0,fill=0.):
        """
        s is the supplied signal part (region)
        start is the time at which it starts
        fill is the value outside the supplied signal
        e.g., s=[1,2,3], start=2, fill=-1
        (starting from time=0)
        [...,-1,-1,1,2,3,-1,...]
        """
        self.s=s
        self.start=start
        self.stop=len(self.s)+start
        self.fill=fill
    def __getitem__(self,i):
        if isinstance(i,slice):
            if i.step and i.step != 1:
                # easiest is just to get the values without step and then do the
                # stepping afterward, but it's not implemented yet
                raise NotImplementedError
            # how long the result is
            len_ret = i.stop - i.start
            # how long the part before region is (can be negative)
            len_pbr = self.start - i.start
            max_0_len_pbr = max(0,len_pbr)
            min_0_len_pbr = min(0,len_pbr)
            # how long the part after the region is (can be negative)
            len_par = i.stop - self.stop
            max_0_len_par = max(0,len_par)
            min_0_len_par = min(0,len_par)
            # how long the inner part is
            len_r = len_ret - (max_0_len_pbr + max_0_len_par)
            # the absolute indices of the inner part
            # start
            start_ir = -min_0_len_pbr
            # stop
            stop_ir = len(self.s) + min_0_len_par
            return np.concatenate((
                np.full(max_0_len_pbr,self.fill,dtype=self.s.dtype),
                self.s[start_ir:stop_ir],
                np.full(max_0_len_par,self.fill,dtype=self.s.dtype),
            ))
        else:
            # assume an integer, lol 
            if i < self.start or i >= self.stop:
                return self.fill
            return self.s[i-self.start]
                
    @property
    def dtype(self):
        return self.s.dtype
