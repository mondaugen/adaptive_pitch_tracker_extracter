# cubic interpolation but vectorized

import numpy as np

def lookup4(x,y):
    """
    x are points to look up
    y are values such that: y[0] = f(0), y[1] = f(1), ..., y[N-1] = f(N-1)
    we must have 1 <= x < N-2
    NOTE: due to the finite precision of x, x will be rounded to the nearest
    floating point number. This can cause intolerable errors in the
    interpolation if x is large. From my experience, for audio signals a good
    rule of thumb is to keep the rounding error less than about 2^-9, if x is
    single precision floats, this corresponds to keeping x within the range x <=
    0 < 2^15 (or 32768). For a larger range a fixed-point implementation can be
    adopted, e.g., Q24.8 would give the same (but constant) worst-case rounding error
    of 2^-9 but have a range of 0 <= x < 2^24.
    """
    N=len(y)
    assert(np.all(1 <= x) and np.all(x < (N-2)))
    x0=np.floor(x).astype('uint32')
    f=x-x0
    f_2=f-2
    f_1=f-1
    f1=f+1
    tmp=np.zeros_like(x)
    ret=np.zeros_like(x)
    ret[:]=y[x0]
    for a in [f_2,f_1,f1,0.5]:
        ret*=a
    tmp[:]=y[x0+1]
    for a in [f_2,f,f1,-0.5]:
        tmp*=a
    ret+=tmp
    tmp[:]=y[x0-1]
    for a in [f,f_1,f_2,-1/6]:
        tmp*=a
    ret+=tmp
    tmp[:]=y[x0+2]
    for a in [f1,f,f_1,1/6]:
        tmp*=a
    ret+=tmp
    return ret
        


