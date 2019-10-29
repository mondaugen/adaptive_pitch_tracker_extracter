# Test avoiding attacks while time stretching using attack_avoider

import numpy as np
import matplotlib.pyplot as plt
import time_map_tstretch
import classic_puckette_timestretch
import common
import attack_finder

SAMPLE_RATE=16000
W=common.get_env('W',default=2048,conv=int)
H=common.get_env('H',default=512,conv=int)
S=common.get_env('S',default=0.5,conv=float)
NG_TH=common.get_env('NG_TH',default=-30,conv=float)

attack_times=np.array([0.5,1.,1.25])
attack_times=(attack_times*SAMPLE_RATE).astype('int')

av=time_map_tstretch.attack_avoider(
attack_times,
-3*H,
W+2*H,
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
x_corr=np.zeros(int(np.round(output_length/S)))
for at,sig in zip(attack_times,sigs):
    x[at:at+len(sig)]+=sig
    ts_at=int(np.round(at/S))
    x_corr[ts_at:ts_at+len(sig)]+=sig
x+=np.random.standard_normal(len(x))*1e-8

analysis_times=np.round(np.arange(0,output_length,H*S)).astype('int')

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

attack_regions=attack_finder.attacks_from_spectral_diff(x,ng_th=NG_TH)
attack_starts=np.array([x for x,y in attack_regions])
attack_ends=np.array([y for x,y in attack_regions])

av_est=time_map_tstretch.attack_avoider(
attack_ends,
-2*H,
W+2*H,
H)
adj_analysis_times_est=[av_est.adjust(at) for at in analysis_times]
reset_times_est = np.where(np.array([b for t,b in adj_analysis_times_est]))[0]
adj_analysis_times_est=np.array([t for t,b in adj_analysis_times_est])
plot_adj_atimes_x_est=np.sort(np.concatenate((adj_analysis_times_est,adj_analysis_times_est)))
ones_osc_est=np.power(-1,np.arange(len(adj_analysis_times_est)))
plot_adj_atimes_y_est=np.zeros(len(plot_adj_atimes_x_est))
plot_adj_atimes_y_est[1::2]=ones_osc_est
plot_adj_atimes_y_est[2::2]=ones_osc_est[:-1]
plot_adj_atimes_y_est[0]=-1

y_est=classic_puckette_timestretch.time_stretch_arb_times(
    x, adj_analysis_times_est, H, W, reset_times=adj_analysis_times_est[reset_times_est])

y_est.tofile('/tmp/attack_avoider_time_stretch_est.f64')

x.tofile('/tmp/attack_avoider_time_stretch_original.f64')
x_corr=x_corr[:len(y)]
x_corr.tofile('/tmp/attack_avoider_time_stretch_correct.f64')
y.tofile('/tmp/attack_avoider_time_stretch_result.f64')
y_no_comp=classic_puckette_timestretch.time_stretch_arb_times(
    x, analysis_times, H, W)
y_no_comp.tofile('/tmp/attack_avoider_time_stretch_no_comp.f64')
print(len(y))
print(len(y_no_comp))

print(reset_times)
fig,axs=plt.subplots(6,1)
axs[0].plot(np.arange(len(x)),x)
axs[0].plot(plot_adj_atimes_x,plot_adj_atimes_y,lw=0.1,c='grey')
axs[0].plot(adj_analysis_times[reset_times],plot_adj_atimes_y[::2][reset_times],'.')
axs[0].set_title('analysis')
axs[1].plot(np.arange(len(y)),y)
axs[1].set_title('synthesis, fixed attacks')
axs[2].plot(np.arange(len(y_no_comp)),y_no_comp)
axs[2].set_title('naive synthesis')
axs[3].plot(np.arange(len(x_corr)),x_corr)
axs[3].set_title('correct synthesis')
axs[4].plot(np.arange(len(x)),x)
print("attack starts")
print(attack_starts)
axs[4].plot(attack_starts,x[attack_starts],'.r')
print("attack ends")
print(attack_ends)
axs[4].plot(attack_ends,x[attack_ends],'.g')
axs[4].plot(plot_adj_atimes_x_est,plot_adj_atimes_y_est,lw=0.1,c='grey')
axs[4].plot(adj_analysis_times_est[reset_times_est],plot_adj_atimes_y_est[::2][reset_times_est],'.')
axs[4].set_title('analysis with estimated attacks')
axs[5].plot(np.arange(len(y_est)),y_est)
axs[5].set_title('synthesis, estimated attacks')

plt.tight_layout()
plt.show()
