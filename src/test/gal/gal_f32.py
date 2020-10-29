# Python bindings for gal_f32.c

from ctypes import *
import numpy as np

lib_gal_f32 = cdll.LoadLibrary("lib_gal_f32.so")

class gal_f32_init(Structure):
    _fields_ = [("P",c_uint)]

class gal_f32_proc(Structure):
    _fields_ = [
    ("x_in",POINTER(c_float)),
    ("_Ef",POINTER(c_float)),
    ("_Eb",POINTER(c_float)),
    ("_R",POINTER(c_float)),
    ("_beta",c_float),
    ("_alpha",c_float),
    ("_D",POINTER(c_float)),
    ("_N",c_uint),
    ("_opt", c_uint),
    ]
    def __init__(self,x,P,opt=0,beta=0.001,alpha=0.001):
        # initialize the pointers using numpy ctypes
        self.x=x.astype('float32')
        self.Ef=np.zeros(len(x),dtype='float32')
        self.Eb=np.zeros(len(x),dtype='float32')
        self.R=np.zeros((len(x),P),dtype='float32')
        self.D=np.ones((len(x),P),dtype='float32')
        super().__init__(
            x_in=self.x.ctypes.data_as(POINTER(c_float)),
            _Ef=self.Ef.ctypes.data_as(POINTER(c_float)),
            _Eb=self.Eb.ctypes.data_as(POINTER(c_float)),
            _R=self.R.ctypes.data_as(POINTER(c_float)),
            _beta=beta,
            _alpha=alpha,
            _D=self.D.ctypes.data_as(POINTER(c_float)),
            _N=len(x),
            _opt=opt,
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

