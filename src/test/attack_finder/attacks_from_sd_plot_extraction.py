import numpy as np
import matplotlib.pyplot as plt

IN_FILE_PATH="/tmp/attacks_from_sd_test_in.f32"
OUT_FILE_PATH_BEG="/tmp/attacks_from_sd_test_out_beg.u32"
OUT_FILE_PATH_END="/tmp/attacks_from_sd_test_out_end.u32"

x=np.fromfile(IN_FILE_PATH,dtype='float32')
x=np.abs(x)
beg=np.fromfile(OUT_FILE_PATH_BEG,dtype='uint32')
print(beg)
end=np.fromfile(OUT_FILE_PATH_END,dtype='uint32')
print(end)

N=len(x)
n=np.arange(N)
plt.plot(n,x)
plt.plot(beg,x[beg],'.',label='beg')
plt.plot(end,x[end],'.',label='end')
plt.legend()

plt.show()
