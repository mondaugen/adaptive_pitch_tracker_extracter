# Time stretch by estimated phase advancement using "ghost" windows.

import numpy as np
from scipy import signal
import common 

def prepare_init_window(w_x_sw,H):
    """
    The first window is the product of the analysis window and the synthesis
    window, overlapped and added so that the initial frames (before the first
    frame) can just be added in without needing to do the whole overlap-and-add
    routine
    """
    W=len(w_x_sw)
    v=np.zeros_like(w_x_sw)
    for h in range(0,W,H):
        v[:W-h] += w_x_sw[h:]
    return v

def tsat_oob_fill_func_default(N):
    return np.random.standard_normal(N)*1e-8

def time_stretch_arb_times(
    x, # signal to time stretch
    t, # times at which to place windows
    H, # hop size
    W, # window size
    window_type='hann',
    synth_window_type=None,
    oob_fill_func=tsat_oob_fill_func_default,
    # if an analysis time is in reset_times, the synthesis frame is taken
    # directly from it and not by summing the phase from the last frame.
    reset_times=[]):

    if synth_window_type is None:
        synth_window_type=window_type

    min_t=np.min(t)
    max_t=np.max(t)

    # pad signal so that all t's are valid
    # beginning of x is considered to be time 0
    t_off=-min(min_t,0)+H
    x=np.concatenate((
        oob_fill_func(t_off).astype(x.dtype),
        x,
        oob_fill_func(max(max_t + W - len(x),0)).astype(x.dtype)))

    # get analysis window
    w=signal.get_window(window_type,W)
    # get synthesis window
    sw=signal.get_window(synth_window_type,W)
    y_len=H*(len(t)-1)+W
    y=np.zeros(y_len)
    # Offset the times by the number of samples padded to the beginning of x
    t += t_off
    # calculate dividing out window, even for hops that aren't multiple of window size
    win_div=1/common.ola_shorten(w*sw,H)
    # calculate initial window so we don't need to do extra ola's to fill the initial buffer
    init_win=prepare_init_window(w*sw,H)
    # sum in "fake" inital overlap-and-add frames
    y[:W]=init_win*x[t[0]:t[0]+W]
    # store initial Y_last
    Y_last=np.fft.fft(x[t[0]:t[0]+W]*w)
    y[:H]*=win_div
    h=H
    for t_ in t[1:]:
        X0=np.fft.fft(x[t_:t_+W]*w)
        if t_-t_off in reset_times:
            Y_last = X0
            y[h:h+W]+=x[t_:t_+W]*w*sw
        else:
            X_H=np.fft.fft(x[t_-H:t_-H+W]*w)
            Y_last=X0*np.abs(X_H)/X_H*Y_last/np.abs(Y_last)
            y[h:h+W]+=np.real(np.fft.ifft(Y_last))*sw
        y[h:h+H]*=win_div
        h+=H

    return y

class ola:
    """
    A port of ola_f32.c
    """
    def _config_chk(sum_in_length,shift_out_length):
        ret = True
        ret &= sum_in_length >= 1
        ret &= shift_out_length <= sum_in_length
        ret &= (common.next_pow_2(shift_out_length + sum_in_length) % shift_out_length) == 0
        if not ret:
            raise ValueError
    def __init__(self,
        # This is usually equal to the window length of the transform
        sum_in_length,
        # This is usually equal to the hop size or equivalently the audio
        # processing block size. This must be a power of 2.
        shift_out_length,
        dtype=np.float64):
        ola._config_chk(sum_in_length,shift_out_length)
        self.sum_in_length = sum_in_length
        self.shift_out_length = shift_out_length
        self.buffer_length = common.next_pow_2(shift_out_length + sum_in_length)
        self.buffer=np.zeros(self.buffer_length,dtype=dtype)
        self.len_mask=self.buffer_length-1
        self.offset = 0
    def sum_in_and_shift_out(self,x_in,overwrite=False):
        """ overwrite forces overwrite, without summing """
        zero_start_0 = (self.offset - self.shift_out_length) & self.len_mask
        if zero_start_0 > self.offset:
            zero_len_0 = self.buffer_length - zero_start_0 
        else:
            zero_len_0 = self.shift_out_length
        zero_len_1 = self.shift_out_length - zero_len_0
        self.buffer[zero_start_0:zero_start_0+zero_len_0] = 0
        self.buffer[:zero_len_1] = 0
        if ((self.offset + self.sum_in_length) & self.len_mask) < self.offset:
            sum_len_0 = self.buffer_length - self.offset
        else:
            sum_len_0 = self.sum_in_length
        sum_len_1 = self.sum_in_length - sum_len_0
        if overwrite:
            self.buffer[self.offset:self.offset+sum_len_0] = x_in[:sum_len_0]
            self.buffer[:sum_len_1] = x_in[sum_len_0:sum_len_0+sum_len_1]
        else:
            self.buffer[self.offset:self.offset+sum_len_0] += x_in[:sum_len_0]
            self.buffer[:sum_len_1] += x_in[sum_len_0:sum_len_0+sum_len_1]
        ret = self.buffer[self.offset:self.offset+self.shift_out_length]
        self.offset = (self.offset + self.shift_out_length) & self.len_mask
        return ret

class pvoc_synth:
    """
    A port of pvoc_synth.c
    """

    def _chk_args(analysis_window,synthesis_window,window_length,hop_size):
        if (window_length < 1) or (hop_size < 1):
            raise ValueError
        if ((len(analysis_window) != window_length) or
                (len(synthesis_window) != window_length)):
            raise ValueError

    def reset_past_output(self):
        """
        Call this to force directly using the input multiplied by the
        init_window. This would be done when you want to simulate a first call
        to process.
        """
        self.z_outputH=None

    def __init__(self,
        analysis_window,
        synthesis_window,
        window_length,
        hop_size,
        # function accepting input_time returning window_length samples
        get_samples):
        pvoc_synth._chk_args(analysis_window,synthesis_window,window_length,hop_size)
        self.analysis_window = analysis_window
        self.synthesis_window = synthesis_window
        self.window_length = window_length
        self.hop_size = hop_size
        self.get_samples = get_samples
        self.ola_buffer = ola(window_length,hop_size)
        self.reset_past_output()
        self.output_scaling = common.ola_shorten(
            self.analysis_window*self.synthesis_window,
            self.hop_size)
        if np.any(self.output_scaling == 0):
            raise ValueError
        self.output_scaling = 1./self.output_scaling
        self.init_window = np.zeros(self.window_length,dtype=self.analysis_window.dtype)
        for h in np.arange(0,self.window_length,self.hop_size):
            a_win_section = self.analysis_window[h:]
            s_win_section = self.synthesis_window[h:]
            self.init_window[:self.window_length-h] += a_win_section*s_win_section

    def process(self,input_time,reset):
        f_input0 = self.get_samples(input_time)
        r_workspace=f_input0*self.analysis_window
        if self.z_outputH is not None:
            z_input0 = np.fft.rfft(r_workspace)
            if reset:
                self.z_outputH = z_input0
                r_workspace *= self.synthesis_window
            else:
                f_inputH = self.get_samples(input_time-self.hop_size)
                r_workspace=f_inputH*self.analysis_window
                z_inputH = np.fft.rfft(r_workspace)
                self.z_outputH *= z_input0*np.abs(z_inputH)/(z_inputH*np.abs(self.z_outputH))
                r_workspace = np.fft.irfft(self.z_outputH)*self.synthesis_window
            self.ola_overwrite=False
        else:
            self.z_outputH=np.fft.rfft(r_workspace)
            r_workspace=f_input0*self.init_window
            self.ola_overwrite=True
        ret = self.ola_buffer.sum_in_and_shift_out(r_workspace,self.ola_overwrite)
        ret *= self.output_scaling
        return ret
