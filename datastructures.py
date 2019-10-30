import numpy as np
import common

class ringbuffer:
    """
    A port of ringbuffer.c
    """
    def __init__(self,length,dtype=np.float64):
        self.buffer_size = common.next_pow_2(length+1)
        self.buffer = np.zeros(self.buffer_size,dtype=dtype)
        self.buf_size_mask = self.buffer_size - 1
        self.head_index = 0
        self.tail_index = 0
    def contents_size(self):
        return (self.tail_idx - self.head_idx) & self.size_mask
    def available_capacity(self):
        return (self.head_idx - 1 - self.tail_idx) & self.size_mask
    def advance_head(self,n):
        """
        Advance the head of the ring buffer. This effectively discards the first n
        values and makes the beginning of the ring buffer n values later than
        previously. If the capacity of the ring buffer is less than n, no advance can be
        made and ValueError is raised.
        """
        if self.contents_size() < n:
            raise ValueError
        self.head_idx = (self.head_idx + n) & self.size_mask
    def push_copy(self,data):
        """
        Put the n values in data at the end of the ring buffer and advance the tail
        index. If there's not enough capacity available, raises ValueError.
        """
        n=len(data)
        if self.available_capacity() < n:
            raise ValueError
        first_region_size = min(self.size - self.tail_idx,n)
        second_region_size = n - first_region_size
        self.buffer[self.tail_idx:self.tail_idx+first_region_size] = data[:first_region_size]
        self.buffer[:second_region_size] = data[first_region_size:]
        self.tail_idx = (self.tail_idx + n) & self.size_mask
    def shift_in(self,data):
        """ Advance the head n values and then push n values in """
        n=len(data)
        self.advance_head(n)
        self.push_copy(data)
    def get_region(self,start,length):
        """
        Start must be in [0,rngbuf_contents_size(rb)-1]
        Length must be in [0,rngbuf_contents_size(rb)-start]
        Makes a copy of the ringbuffer contents and returns them.
        """
        contents_size = self.contents_size()
        if start >= contents_size:
            raise IndexError
        if length > (contents_size-start):
            raise ValueError
        start_idx = (start + self.head_idx) & self.size_mask
        ret=np.zeros(length,dtype=self.buffer.dtype)
        first_region_size = min(self.size - start_idx,length)
        second_region_size = length - first_region_size
        ret[:first_region_size]=self.data[start_idx:start_idx+first_region_size]
        ret[first_region_size:]=self.data[:second_region_size]
        return ret
