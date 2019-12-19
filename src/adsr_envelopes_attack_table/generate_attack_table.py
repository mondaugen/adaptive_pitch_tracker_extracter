import numpy as np
from scipy import signal
import subprocess
import os

def is_pow_2(x):
    t=1
    while (t < x):
        t <<= 1
    return t == x

file_stem=os.environ['OUTPUT_FILE_STEM']
# this has to be a power of 2
N=4096
assert(is_pow_2(N))
table=signal.get_window('hann',2*N)[:N]
table=table.astype('float32')
# get dirname $0
dn=os.path.dirname(__file__)
table.tofile(os.path.join(dn,file_stem+'.f32'))
with open(os.path.join(dn,file_stem+'.h'),'w') as f:
    f.write("#define %s_length %d\n" % (file_stem,N,))
    f.write("extern const float _binary_%s_f32_start;\n" % (file_stem,))
    f.write("#define %s (&_binary_%s_f32_start)\n" % (file_stem,file_stem,))
    

