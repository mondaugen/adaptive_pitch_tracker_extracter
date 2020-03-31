import numpy as np
import matplotlib.pyplot as plt
import common

IN_FILE=common.get_env('IN_FILE',default='/tmp/in.f64')
ATTACKS_FILE=common.get_env('ATTACKS_FILE',default='/tmp/attacks.u32')
LOOK_UP_TIMES=common.get_env('LOOK_UP_TIMES',default='/tmp/look_up_times.u32')
RESET_TIMES=common.get_env('RESET_TIMES',default='/tmp/reset_times.u32')
SR=common.get_env('SR',default=16000,conv=float)
H=common.get_env('H',default=1024,conv=int)
W=common.get_env('W',default=4096,conv=int)

x=np.fromfile(IN_FILE,dtype='float64')
attacks=np.fromfile(ATTACKS_FILE,dtype='uint32')
look_up_times=np.fromfile(LOOK_UP_TIMES,dtype='uint32')
reset_times=np.fromfile(RESET_TIMES,dtype='uint32')

# append 0s to x because it is possible we made a look up outside of its valid
# domain
look_up_times=look_up_times[look_up_times>=0]
reset_times=reset_times[reset_times>=0]
pad=np.zeros(max(0,1+max(max(look_up_times),max(reset_times))-len(x)),dtype=x.dtype)
x=np.concatenate((x,pad))

N=len(x)
n=np.arange(len(x))
t=n/SR
fig,ax=plt.subplots(1,1)
ax.plot(t,x)
ax.plot(t[attacks],x[attacks],'.',label='attacks')
z=np.zeros_like(x)
z[look_up_times]=1
z[reset_times]+=1
resets=(z[look_up_times]>1)
h=1
for l,r in zip(look_up_times,resets):
    common.plot_arch(ax,l/SR,W/SR,h,color='red' if r else 'black',print_h='%d',t=W/SR/16)
    h+=1

plt.show()
