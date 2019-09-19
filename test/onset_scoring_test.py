# Use mir_eval to evaluate the onset score

import synth_af_score
import attack_finder
import numpy as np
import mir_eval
import common
import matplotlib.pyplot as plt

def gef(n,f):
    return common.get_env(n,default=f,conv=float)

def gei(n,i):
    return common.get_env(n,default=i,conv=int)

SOUNDFILE=common.get_env('SOUNDFILE',check_if_none=True)
SAMPLE_RATE=common.get_env('SAMPLE_RATE',default=16e3,conv=float)
SCOREFILE=common.get_env('SCOREFILE',check_if_none=True)

# true time stamps in seconds
ts=synth_af_score.extract_ts(SCOREFILE)
ts=np.sort(np.array(ts))

x=np.fromfile(SOUNDFILE)
# estimated time stamps
H=gei('H',256)
W=gei('W',256)
alpha=gef('ALPHA',1)
thresh=gef('THRESH',3e-3)
ts_est,a=attack_finder.attack_region_estimation(x,H,W,alpha,thresh)
a_ts=ts_est[:-1][(np.diff(a)>0)]/SAMPLE_RATE


XLIM=gef('XLIM',2)

fig,axs=plt.subplots(3,1)
t=np.arange(len(x))/SAMPLE_RATE
axs[0].plot(t,x)
axs[0].set_title('Original')
axs[1].plot(ts_est/SAMPLE_RATE,a)
axs[1].set_title('Attack regions')
axs[2].plot(ts_est[:-1]/SAMPLE_RATE,np.array(np.diff(a)>0,dtype='float'))
axs[2].set_title('Attack points')

for ax in axs:
    ax.set_xlim(0,XLIM)

print("True")
print(ts)
print("Est")
print(a_ts)

scores=mir_eval.onset.evaluate(ts,a_ts)
for key in scores.keys():
    print(key,scores[key])

plt.show()
