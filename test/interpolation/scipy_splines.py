import common
import numpy as np
from scipy.interpolate import interp1d

INPUT=common.get_env('INPUT',default='/tmp/in.f64')
OUTPUT=common.get_env('OUTPUT',default='/tmp/out.f64')
S=common.get_env('S',default=1.,conv=float)

x=np.fromfile(INPUT,dtype='float64')
n=np.arange(len(x))
t=np.arange(0,len(x)-1,S)

interp=interp1d(n,x,kind='cubic')
y=interp(t)

y.tofile(OUTPUT)

