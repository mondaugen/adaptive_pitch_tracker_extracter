import numpy as np
from wavetables import fast_log2_aprox_frac_f32

print(2,fast_log2_aprox_frac_f32(np.array(2,dtype='float32')),np.log2(2))
print(3,fast_log2_aprox_frac_f32(np.array(3,dtype='float32')),np.log2(3))
print(4,fast_log2_aprox_frac_f32(np.array(4,dtype='float32')),np.log2(4))
print(5,fast_log2_aprox_frac_f32(np.array(5,dtype='float32')),np.log2(5))
print(6,fast_log2_aprox_frac_f32(np.array(6,dtype='float32')),np.log2(6))
print(7,fast_log2_aprox_frac_f32(np.array(7,dtype='float32')),np.log2(7))
print(8,fast_log2_aprox_frac_f32(np.array(8,dtype='float32')),np.log2(8))
print(12,fast_log2_aprox_frac_f32(np.array(12,dtype='float32')),np.log2(12))
