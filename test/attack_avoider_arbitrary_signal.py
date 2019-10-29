# Test avoiding attacks while time stretching using attack_avoider

import numpy as np
import matplotlib.pyplot as plt
import time_map_tstretch
import classic_puckette_timestretch
import common


SAMPLE_RATE=16000
W=2048
H=256

attack_times=np.array([0.5,1.,1.5])
attack_times=(attack_times*SAMPLE_RATE).astype('int')

av=time_map_tstretch.attack_avoider(
attack_times,
-4*H,
W+3*H,
H)

# make signal

sigs=[]
for fname in [
    '/tmp/piano_adj_beg/60.f64',
    '/tmp/piano_adj_beg/64.f64',
    '/tmp/piano_adj_beg/67.f64']:
    sigs.append(np.fromfile(fname))

file_lengths=np.array([len(s) for s in sigs])
output_length=np.max(np.add.outer(attack_times,file_lengths))
x=np.zeros(output_length)
for at,sig in zip(attack_times,sigs):
    x[at:at+len(sig)]+=sig
x+=np.random.standard_normal(len(x))*1e-8

# generate some interesting analysis times
n_rates=100
rates=np.random.standard_normal(n_rates)*2
interp_rates=np.interp(
np.arange(output_length),
np.linspace(0,output_length,n_rates),
rates)
interp_pos=np.cumsum(interp_rates)
analysis_times=np.round(np.mean(common.frame(interp_pos,H,H),axis=0)).astype('int')

adj_analysis_times=[av.adjust(at) for at in analysis_times]
reset_times = np.where(np.array([b for t,b in adj_analysis_times]))[0]
adj_analysis_times=np.array([t for t,b in adj_analysis_times])
plot_adj_atimes_x=np.sort(np.concatenate((adj_analysis_times,adj_analysis_times)))
ones_osc=np.power(-1,np.arange(len(adj_analysis_times)))
plot_adj_atimes_y=np.zeros(len(plot_adj_atimes_x))
plot_adj_atimes_y[1::2]=ones_osc
plot_adj_atimes_y[2::2]=ones_osc[:-1]
plot_adj_atimes_y[0]=-1

print(plot_adj_atimes_x[:10])
print(plot_adj_atimes_y[:10])


y=classic_puckette_timestretch.time_stretch_arb_times(
    x, adj_analysis_times, H, W, reset_times=adj_analysis_times[reset_times])

x.tofile('/tmp/attack_avoider_arb_sig_original.f64')
y.tofile('/tmp/attack_avoider_arb_sig_result.f64')
y_no_comp=classic_puckette_timestretch.time_stretch_arb_times(
    x, analysis_times, H, W)
y_no_comp.tofile('/tmp/attack_avoider_arb_sig_no_comp.f64')
print(len(y))
print(len(y_no_comp))


print(reset_times)
fig,axs=plt.subplots(4,1)
axs[0].plot(np.arange(len(x)),x)
#axs[0].plot(plot_adj_atimes_x,plot_adj_atimes_y)
#axs[0].plot(adj_analysis_times[reset_times],plot_adj_atimes_y[::2][reset_times],'.')
axs[0].set_title('analysis')
axs[1].plot(np.arange(len(y)),y)
axs[1].set_title('synthesis, fixed attacks')
axs[2].plot(np.arange(len(y_no_comp)),y_no_comp)
axs[2].set_title('naive synthesis')
axs[3].plot(np.arange(output_length),interp_rates)

plt.tight_layout()
plt.show()
