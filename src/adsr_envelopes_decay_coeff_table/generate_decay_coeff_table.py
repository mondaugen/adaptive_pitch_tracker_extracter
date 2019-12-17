import numpy as np
from scipy import signal
import subprocess
import os

file_stem=os.environ['OUTPUT_FILE_STEM']
min_N=48000*5
decay_min_dB=-90
decay_min_A=np.power(10,decay_min_dB/20)
tab_len=0
# We give the table an extra point at the end because the
# adsr_decay_coeff_lookup function restricts the argument in [1,2**(tab_len-2)]
while (1 << max(tab_len-1,0)) < min_N:
    tab_len+=1
table_x=np.power(2,np.arange(tab_len))
table=np.power(decay_min_A,1/table_x)
table=table.astype('float32')
# get dirname $0
dn=os.path.dirname(__file__)
table.tofile(os.path.join(dn,file_stem+'.f32'))
with open(os.path.join(dn,file_stem+'.h'),'w') as f:
    f.write("#define %s_length %d\n" % (file_stem,tab_len,))
    f.write("#define adsr_decay_min_A %.18g\n" % (decay_min_A,))
    f.write("extern const float _binary_%s_f32_start;\n" % (file_stem,))
    f.write("#define %s (&_binary_%s_f32_start)\n" % (file_stem,file_stem,))
    

