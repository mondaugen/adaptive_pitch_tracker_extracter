import datastructures as ds
import numpy as np

rb1=ds.ringbuffer(10)

assert(rb1.contents_size() == 0)

assert(rb1.available_capacity() >= 10)

flg=0
try:
    rb1.advance_head(5)
except ValueError:
    flg=1
assert(flg)

rb1.push_copy(np.arange(10))
    
rb1.shift_in(np.arange(5))

x=rb1.get_region(0,10)
assert(np.all(x==np.concatenate((np.arange(10)[5:],np.arange(5)))))
print(x)

flg=0
rb1.advance_head(5)
try:
    x=rb1.get_region(0,10)
except ValueError:
    flg=1
assert(flg)

flg=0
try:
    x=rb1.get_region(5,1)
except IndexError:
    flg=1
assert(flg)

x=rb1.get_region(0,5)
assert(np.all(x==np.arange(5)))
print(x)
