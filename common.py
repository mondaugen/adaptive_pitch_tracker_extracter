import numpy as np

def normalize(x):
    x-=np.mean(x)
    x/=np.max(np.abs(x))
    return x

