from ctypes import *
import numpy as np

lib_iirl_f32 = cdll.LoadLibrary("lib_iirl_f32.so")

class iirlf_f32_init(Structure):
    _fields_ = [("P",c_uint)]

class iirlf_f32_proc(Structure):
    _fields_ = [("xin",POINTER(c_float)),
                ("out",POINTER(c_float)),
                ("N",c_uint),
                ("R",POINTER(c_float)),
                ("b0",c_float),
                ("opts",c_uint)]
    def __init__(self,xin,R,b0,varying_R=True):
        self._xin=xin.astype("float32")
        self._R=R.astype("float32")
        self._b0=b0
        self._out=np.zeros_like(self._xin)
        self._opts = 1 if varying_R else 0
        super().__init__(
            xin=self._xin.ctypes.data_as(POINTER(c_float)),
            out=self._out.ctypes.data_as(POINTER(c_float)),
            N=len(self._xin),
            R=self._R.ctypes.data_as(POINTER(c_float)),
            b0=self._b0,
            opts=self._opts
        )

class iirlf_f32:
    def __init__(self,P):
        init = iirlf_f32_init(P)
        self.handle = lib_iirl_f32.iir_lattice_filter_f32_new(byref(init))
        if self.handle == c_voidp(0):
            raise ValueError("Couldn't initialize with P=%d" %(P,))
    def __del__(self):
        if self.handle:
            lib_iirl_f32.iir_lattice_filter_f32_free(self.handle)
    def proc(self,s):
        lib_iirl_f32.iir_lattice_filter_f32_process(self.handle, byref(s))

