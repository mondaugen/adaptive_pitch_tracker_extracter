import numpy as np
#from scipy import signal
import subprocess
import os

# compilation will be fast even with N so large
N=4000000
#table=signal.get_window('hann',2*N)[1:N+1]
table=np.arange(N)
table=table.astype('float32')
# get dirname $0
dn=os.path.dirname(__file__)
table.tofile(os.path.join(dn,'table.f32'))
with open(os.path.join(dn,'table.h'),'w') as f:
    f.write("#define table_N %d\n" % (N,))
    f.write("extern const float _binary___table_f32_start;\n")
    f.write("#define table (&_binary___table_f32_start)\n")
    
