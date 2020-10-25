import filters
import gal_alexander
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import os
import common

P=int(os.environ.get("P","50"))
INFILE=os.environ.get("INFILE","/tmp/in.f32")
OUTFILE=os.environ.get("OUTFILE","/tmp/out.f32")
OUTFILE_FILTERED=os.environ.get("OUTFILE_FILTERED","/tmp/out_filt.f32")
x=np.fromfile(INFILE,dtype=np.float32)
N=len(x)
Ef,Eb,D,K=gal_alexander.ngal(x,P,alpha=0.001,beta=0.001,normalize=True)
Ef[1:N+1,-1].astype('float32').tofile(OUTFILE)
zz=np.random.standard_normal(N)
yy=filters.iir_lattice_filter_proc(zz,K[1:N+1,:].T)
yy=common.normalize(yy)
yy.astype('float32').tofile(OUTFILE_FILTERED)
