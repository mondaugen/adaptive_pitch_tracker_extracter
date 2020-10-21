# Python bindings for gal_f32.c

from ctypes import *
import numpy as np

lib_gal_f32 = cdll.LoadLibrary("lib_gal_f32.so")

class gal_f32_init(Structure):
    _fields_ = [("P",c_uint)]

class gal_f32_proc(Structure):
    _fields_ = [
    ("x_in",POINTER(c_float)),
    ("_ep",POINTER(c_float)),
    ("_em",POINTER(c_float)),
    ("_R",POINTER(c_float)),
    ("_mu",POINTER(c_float)),
    ("_beta",c_float),
    ("_l",c_float),
    #TODO: Can we properly represent an enum?
    ("_opt", c_uint),
    ("_N",c_uint)
    ]
    def __init__(self,x,mu):
        # initialize the pointers using numpy ctypes
        self.x=x.astype('float32')
        self.ep=np.zeros(len(x),dtype='float32')
        self.em=np.zeros(len(x),dtype='float32')
        self.R=np.zeros((len(x),len(mu)),dtype='float32')
        self.mu=mu.astype('float32')
        super().__init__(
            x_in=self.x.ctypes.data_as(POINTER(c_float)),
            _ep=self.ep.ctypes.data_as(POINTER(c_float)),
            _em=self.em.ctypes.data_as(POINTER(c_float)),
            _R=self.R.ctypes.data_as(POINTER(c_float)),
            _mu=self.mu.ctypes.data_as(POINTER(c_float)),
            _beta=0,
            _l=0,
            _opt=0,
            _N=len(x)
        )



class gal_f32:
    def __init__(self,P):
        init = gal_f32_init(P)
        self.handle = lib_gal_f32.gal_f32_new(byref(init));
        if self.handle == c_voidp(0):
            raise ValueError("Couldn't initialize with P=%d" %(P,))
    def __del__(self):
        if self.handle:
            lib_gal_f32.gal_f32_free(self.handle)
    def proc(self,s):
        lib_gal_f32.gal_f32_process(self.handle, byref(s))

