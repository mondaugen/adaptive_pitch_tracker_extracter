from ctypes import *
import numpy as np
import sys

libadder = CDLL("libadder.so")
print(libadder.adder_add)

class adder_proc:
    class _ctype(Structure):
        _fields_ = [
            ("N",c_uint),
            ("a",POINTER(c_float)),
            ("b",POINTER(c_float))
        ]
    def __init__(self,N,a,b):
        self.a=np.array(a,dtype='float32')
        self.b=np.array(b,dtype='float32')
        self.struct = adder_proc._ctype(N,
            self.a.ctypes.data_as(POINTER(c_float)),
            self.b.ctypes.data_as(POINTER(c_float)))

ap=adder_proc(
5,
[1,2,3,4,5],
[5,4,3,2,1]
)

libadder.adder_add(byref(ap.struct))

print(ap.a)
print(ap.b)
