import numpy as np
import common

class ringbuffer:
    """
    A port of ringbuffer.c
    """
    def reset(self):
        """ Set head_index and tail_index to 0 making the ringbuffer empty. """
        self.head_index = 0
        self.tail_index = 0
    def __init__(self,length,dtype=np.float64):
        self.buffer_size = common.next_pow_2(length+1)
        self.buffer = np.zeros(self.buffer_size,dtype=dtype)
        self.buf_size_mask = self.buffer_size - 1
        self.reset()
    def contents_size(self):
        return (self.tail_index - self.head_index) & self.buf_size_mask
    def available_capacity(self):
        return (self.head_index - 1 - self.tail_index) & self.buf_size_mask
    def advance_head(self,n):
        """
        Advance the head of the ring buffer. This effectively discards the first n
        values and makes the beginning of the ring buffer n values later than
        previously. If the capacity of the ring buffer is less than n, no advance can be
        made and ValueError is raised.
        """
        if self.contents_size() < n:
            raise ValueError
        self.head_index = (self.head_index + n) & self.buf_size_mask
    def push_copy(self,buffer):
        """
        Put the n values in buffer at the end of the ring buffer and advance the tail
        index. If there's not enough capacity available, raises ValueError.
        """
        n=len(buffer)
        if self.available_capacity() < n:
            raise ValueError
        first_region_size = min(self.buffer_size - self.tail_index,n)
        second_region_size = n - first_region_size
        self.buffer[self.tail_index:
            self.tail_index+first_region_size] = buffer[:first_region_size]
        self.buffer[:second_region_size] = buffer[first_region_size:]
        self.tail_index = (self.tail_index + n) & self.buf_size_mask
    def shift_in(self,buffer):
        """ Advance the head n values and then push n values in """
        n=len(buffer)
        self.advance_head(n)
        self.push_copy(buffer)
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
        start_index = (start + self.head_index) & self.buf_size_mask
        ret=np.zeros(length,dtype=self.buffer.dtype)
        first_region_size = min(self.buffer_size - start_index,length)
        second_region_size = length - first_region_size
        ret[:first_region_size]=self.buffer[start_index:start_index+first_region_size]
        ret[first_region_size:]=self.buffer[:second_region_size]
        return ret
