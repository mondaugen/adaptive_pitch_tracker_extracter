import numpy as np

x=np.arange(-100000,100000,dtype='int64')<<16
x.tofile("/tmp/dspm_2dline_s48q16_lookup_vs48q16_input.s48q16")
