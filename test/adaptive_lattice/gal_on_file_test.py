import filters
import gal_f32
from iir_lattice_filter import iirlf_f32, iirlf_f32_proc
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import os
import common

P=int(os.environ.get("P","300"))
INFILE=os.environ.get("INFILE","/tmp/in.f32")
OUTFILE=os.environ.get("OUTFILE","/tmp/out.f32")
OUTFILE_FILTERED=os.environ.get("OUTFILE_FILTERED","/tmp/out_filt.f32")
OUTFILE_NOISE_FILTERED=os.environ.get("OUTFILE_NOISE_FILTERED","/tmp/out_noise_filt.f32")
x=np.fromfile(INFILE,dtype=np.float32)
N=len(x)
gal_proc=gal_f32.gal_f32_proc(x,P,alpha=0.3,beta=0.3,opt=0)
gal=gal_f32.gal_f32(P)
gal.proc(gal_proc)
gal_proc.Ef.tofile(OUTFILE)
iirlf=iirlf_f32(P)
# We can get very close but not when normalize=True
iirlfp=iirlf_f32_proc(gal_proc.Ef,gal_proc.R,1)
iirlf.proc(iirlfp)
print(np.mean(np.abs(iirlfp._out-x)))
yy=iirlfp._out
yy=common.normalize(yy)
yy.astype('float32').tofile(OUTFILE_FILTERED)

zz=np.random.standard_normal(N)
iirlf=iirlf_f32(P)
iirlfp=iirlf_f32_proc(zz,gal_proc.R,1)
iirlf.proc(iirlfp)
yy=iirlfp._out
yy=common.normalize(yy)
yy.astype('float32').tofile(OUTFILE_NOISE_FILTERED)
