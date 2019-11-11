import numpy as np
from datastructures import ringbuffer

# Warnings
# When time is adjusted to be in bounds
TIME_ADJUSTED=(1 << 0)
# When length is adjusted to not exceed delay line length
LEN_ADJUSTED=(1 << 1)

class access_struct:
    """ This gets subclassed and the subclass implements the __call__ method.
    The subclass can use the fields of an instance of this class (they will be
    initialized by routines). The caller can also add additional fields to this
    subclass if it needs them but note that callers will generally overwrite the
    superclass rel_del_line_access's fields. """
    def __init__(self,values=None,warnings=0):
        self.values = values
        self.warnings = warnings

class rel_del_line:
    """ A delay line whose values can be looked up relative to a time. """
    def reset(self):
        self.rb.push_copy(np.zeros(self.length,dtype=self.dtype))
    def __init__(self,length,dtype=np.float64,print_warnings=False):
        self.dtype = dtype
        self.time_zero_ago = 0
        self.rb = ringbuffer(length,dtype=dtype)
        self.length = length
        self.print_warnings = print_warnings
        self.reset()
    def access(self,
        # The largest valid time is (time + length - 1) <= (self.time_zero_ago - 1) or
        # time <= self.time_zero_ago - length
        # The smallest valid time is self.time_zero_ago - self.length
        # See below for what happens if these rules are not respected.
        time,
        length,
        # subclass of access_struct at least implementing __call__ method
        access_fun):
        access_fun.warnings = 0
        if time < (self.time_zero_ago - self.length):
            # adjust time, give warning
            time = self.time_zero_ago  - self.length + 1
            access_fun.warnings |= TIME_ADJUSTED
            if self.print_warnings:
                print("TIME_ADJUSTED")
        if time > (self.time_zero_ago - length):
            # adjust length, give warning
            length = self.time_zero_ago - time
            access_fun.warnings |= LEN_ADJUSTED
            if self.print_warnings:
                print("LEN_ADJUSTED")
        look_up_time = time - self.time_zero_ago
        region = self.rb.get_region(self.rb.contents_size() + look_up_time, length)
        access_fun.values = region
        access_fun()
    def process(self, values):
        len_values=len(values)
        self.rb.shift_in(values)
        self.time_zero_ago += len_values
