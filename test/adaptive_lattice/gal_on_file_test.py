import filters
import gal_f32
from iir_lattice_filter import iirlf_f32, iirlf_f32_proc
import numpy as np
from scipy import signal, interpolate
import matplotlib.pyplot as plt
import os
import common

P=int(os.environ.get("P","300"))
# time stretch amount
S=float(os.environ.get("S","1"))
INFILE=os.environ.get("INFILE","/tmp/in.f32")
OUTFILE=os.environ.get("OUTFILE","/tmp/out.f32")
OUTFILE_FILTERED=os.environ.get("OUTFILE_FILTERED","/tmp/out_filt.f32")
OUTFILE_NOISE_FILTERED=os.environ.get("OUTFILE_NOISE_FILTERED","/tmp/out_noise_filt.f32")
OUTFILE_TS=os.environ.get("OUTFILE_TS","/tmp/out_ts.f32")
OUTFILE_FILTERED_TS=os.environ.get("OUTFILE_FILTERED_TS","/tmp/out_filt_ts.f32")
x=np.fromfile(INFILE,dtype=np.float32)
N=len(x)
n=np.arange(N)
gal_proc=gal_f32.gal_f32_proc(x,P,alpha=0.01,beta=0.1,opt=0)
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

# time stretch
nS=np.arange(0,N-1,S)
# time stretched signal
Ef_ts=interpolate.interp1d(n,gal_proc.Ef)(nS).astype('float32')
Ef_ts.tofile(OUTFILE_TS)
# time stretched R
R_ts=interpolate.interp1d(n,gal_proc.R,axis=0)(nS).astype('float32')
QQ=np.zeros_like(Ef_ts)
QQ[np.random.randint(0,np.floor(N/S)-1,100)]=1
iirlf=iirlf_f32(P)
#iirlfp=iirlf_f32_proc(Ef_ts,R_ts,1)
iirlfp=iirlf_f32_proc(QQ,R_ts,1)
iirlf.proc(iirlfp)
xx=common.normalize(iirlfp._out).astype('float32')
xx.tofile(OUTFILE_FILTERED_TS)
