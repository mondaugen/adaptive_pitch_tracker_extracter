import numpy as np
import matplotlib.pyplot as plt
import spectral_difference
import os

env=os.environ

one_sided_max="both"
try:
    one_sided_max=env['ONE_SIDED_MAX']
except KeyError:
    pass

truth_src=np.fromfile("/tmp/local_max_f32_input.f32",dtype='float32')
truth=spectral_difference.local_max(truth_src,one_sided_max=one_sided_max)
test=np.fromfile( "/tmp/local_max_f32_output.f32",dtype='uint32')
try:
    _ = env['NO_TEST']
except KeyError:
    if np.all((truth-test)==0):
        print("Passed")
    else:
        print("Failed")
N=len(truth_src)
n=np.arange(N)
plt.plot(n,truth_src)
plt.plot(truth,truth_src[truth],'o',label='truth')
plt.plot(test,truth_src[test],'.',label='test')
plt.legend()
plt.show()
