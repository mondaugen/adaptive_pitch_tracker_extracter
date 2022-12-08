from freq_dom_window import freq_dom_window, dft_dv, calc_X, dft, sum_of_cos_dft_win_type
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from dftdk import gradient_ascent_step_harm_lock
from some_ft import normalize_sum_of_cos_A, sum_of_cos_dft
from some_sig import comb_no_mod, sum_of_cos
from scipy import signal, interpolate
from fdw_tracker import fdw_tracker, v_dev_stop_criterion, t_bounds_stop_criterion, multi_analyse_combine
from cubic_sinusoid_synth import quadratic_phase_poly_interp, linear_amplitude_interp
from common import normalize, get_env
import partial_processing as pp
from buf import inf_buf

IN_FILE=get_env('IN_FILE',default='sounds/tanguillo-1ch-48k.f64')
OUT_FILE=get_env('OUT_FILE',default='/tmp/a.f64')
SR=get_env('SR',default=44100,conv=float)
T_START=get_env('T_START',default=0.6,conv=float)
F_START=get_env('F_START',default=150.,conv=float)
R=get_env('R',default=0.02,conv=float) # output rate
# amount of time at the beginning that no time stretching is performed
T_NO_R=get_env('T_NO_R',conv=float)
# stop criterion
SC=get_env('SC',default='v_dev')
# v_dev stop criterion parameters
MAX_DEV_SEMITONES=get_env('MAX_DEV_SEMITONES',default=0.5,conv=float)
# t_bounds stop criterion parameters
T_MIN=get_env('T_MIN',conv=float)
T_MAX=get_env('T_MAX',conv=float)


x=np.fromfile(IN_FILE)
sr=SR
t_start=T_START
n_start=int(np.round(t_start*sr))

# starting frequency
vstart=np.array([F_START/sr])
n_harms=30
vstarts=np.multiply.outer(vstart,(1+np.arange(n_harms))).flatten()
v_groups=np.multiply.outer(np.arange(len(vstart)),np.ones(n_harms)).flatten().astype('int')
print(v_groups)
N=2048
N_h=N//8
fdwt=fdw_tracker(N=N)

def sc():
    if SC == 'v_dev':
        return v_dev_stop_criterion(vstarts[0],max_dev_semitones=MAX_DEV_SEMITONES)
    if SC == 't_bounds':
        h_min=int(round(T_MIN*sr))
        h_max=int(round(T_MAX*sr))
        return t_bounds_stop_criterion(h_min,h_max)

wrapped_x=inf_buf(x)
h_fw=np.arange(n_start,len(x),N_h)
h_bw=np.arange(n_start,0,-N_h)

k,h,X,gr,X_fr,k_fr = multi_analyse_combine(fdwt,x=wrapped_x,vstarts=vstarts,h=[h_bw,h_fw],v_groups=v_groups,warm_up_hops=0,mu=0.2,n_steps=3,grad_weight='equal',stop_criterion=sc)
gr_sc=np.exp(-10*(gr*gr))
fig_spec,ax_spec=plt.subplots()
def plot_spec_line(line,col):
    x=k_fr[:,col]/N*SR
    y=20*np.log10(np.abs(X_fr[:,col]))
    if line is None:
        line, = ax_spec.plot(x,y)
    else:
        line.set_ydata(y)
    return line
def plot_ana_points(line,col):
    x=k[:,col]/N*SR
    y=20*np.log10(np.abs(X[:,col]))
    if line is None:
        line, = ax_spec.plot(x,y,'.')
    else:
        line.set_xdata(x)
        line.set_ydata(y)
    return line
def plot_amp_fit_points(line,col):
    def regX(x):
        return np.hstack((np.ones(len(x))[:,None],x[:,None],x[:,None]**2))
    x=k[:,col]/N*SR
    A=regX(x)
    B=20*np.log10(np.abs(X[:,col]))
    W=np.diag(gr_sc[:,col])
    prms=linalg.lstsq(W@A,W@B)
    x_plt=np.arange(0,N,1./8)/N*SR
    y_plt=regX(x_plt)@prms[0]
    if line is None:
        line, = ax_spec.plot(x_plt,y_plt)
    else:
        line.set_xdata(x_plt)
        line.set_ydata(y_plt)
    return line

    
    

ax_spec_slider = fig_spec.add_axes([0.25, 0.1, 0.65, 0.03])
spec_slider=Slider(ax_spec_slider,'Frame',0,X_fr.shape[1]-1,valinit=0,valstep=1)
class update_spec:
    def __init__(self):
        self.spec_line = plot_spec_line(None,0)
        self.ana_points = plot_ana_points(None,0)
        self.amp_fit_points = plot_amp_fit_points(None,0)
        ax_spec.set_ylim([-100,0])
    def __call__(self,val):
        self.spec_line = plot_spec_line(self.spec_line,spec_slider.val)
        self.ana_points = plot_ana_points(self.ana_points,spec_slider.val)
        self.amp_fit_points = plot_amp_fit_points(self.amp_fit_points,spec_slider.val)
        ax_spec.set_ylim([-100,0])
        fig_spec.canvas.draw_idle()
spec_slider.on_changed(update_spec())
    
h_plot = h

plt.figure(2)
plot_tracks=True
if plot_tracks:
    for k_row in k:
        plt.plot(h_plot/sr,k_row[:len(h)]/N*sr)

plt.specgram(x,NFFT=N,Fs=sr,
window=fdwt.win_type.window(N,centered_at_time_zero=False),
noverlap=N-np.abs(N_h),sides='onesided')

X_synth=X
X_abs=np.abs(X_synth)
X_phs=np.angle(X_synth)
omega=k/N*2*np.pi
if N_h < 0:
    omega=omega[:,::-1]

# output rate for time stretching (<1) or compression (>1)
output_rate=R
if T_NO_R is None:
    orig_rate_len=omega.shape[1]
else:
    nh_no_r=int(round(T_NO_R*sr/N_h))
    orig_rate_len=min(nh_no_r,omega.shape[1])
input_indices=np.arange(0,omega.shape[1])
output_indices=np.concatenate((
    np.arange(0,orig_rate_len),
    np.arange(orig_rate_len,omega.shape[1]-1+output_rate,output_rate))
)
# It seems sometimes output_indices can exceed the valid boundaries, maybe
# because of accumulative errors in np.arange
output_indices=output_indices[(0 <= output_indices) & (output_indices <= (omega.shape[1]-1))]
#X_abs=pp.squish_anomalies(X_abs)

phase_track_interp=interpolate.interp1d(
input_indices,
omega,
axis=1)
omega_interp=phase_track_interp(output_indices)

amp_track_interp=interpolate.interp1d(
input_indices,
X_abs,
axis=1)
X_abs_interp=amp_track_interp(output_indices)

phase_tracks=quadratic_phase_poly_interp(np.abs(N_h),X_phs[:,0],omega_interp.T)
amp_tracks=linear_amplitude_interp(np.abs(N_h),X_abs_interp.T)
print('amp_tracks.shape',amp_tracks.shape)
print('X_abs.shape',X_abs.shape)
fig,ax=plt.subplots()
ax_track = fig.add_axes([0.25, 0.05, 0.65, 0.03])
ax_gr=ax.twinx()
# the figure 2 plotting function
t_amp_track=h_plot/sr
def track_line(p):
    track = 20*np.log10(X_abs[p,:])
    return track
tr_line,=ax.plot(t_amp_track,track_line(0))
gr_line,=ax_gr.plot(t_amp_track,gr_sc[0,:],c='g')
ax.set_ylim([-100,0])
ax_gr.set_ylim([0,1])
track_slider=Slider(ax_track,"Track",0,X_abs.shape[0]-1,valinit=0,valstep=1)
def update_fig2(val):
    tr_line.set_ydata(track_line(track_slider.val))
    gr_line.set_ydata(gr_sc[track_slider.val,:])
    ax.set_ylim([-100,0])
    fig.canvas.draw_idle()
track_slider.on_changed(update_fig2)

x_synth=np.real((np.exp(complex('j')*phase_tracks)*amp_tracks).sum(axis=0))
normalize(x_synth).tofile(OUT_FILE)

plt.show()
