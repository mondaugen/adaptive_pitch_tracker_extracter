# A pitch shifter that works by obtaining a signal from something (e.g., a time
# stretcher) and resamples this signal.
from datastructures import ringbuffer
from scipy import interpolate
import numpy as np

def default_get_interpolator(x,y):
    # Note this interpolator is not exactly what we want to use in practice. For
    # signal indices in the range [0,...,N-1], we will use N Lagrange
    # polynomials computed using the values [-1,0,1,2], [0,1,2,3], ...,
    # [N-2,N-1,N,N+1]. That means we need -1, N and N+1. In other words, if x[0]
    # represents index -1 etc, then this function should fail if you input a
    # value between x[0] and x[1], but that is not the case for the scipy
    # implementation.
    return interpolate.interp1d(x,y,kind='cubic',assume_sorted=True)

def default_get_interpolator_range(x_min,x_max):
    return (int(np.floor(x_min)-1),int(np.floor(x_max)+2))

def default_get_interpolator_n_points(N):
    return N+3

class pitch_shifter:
    def __init__(self,
        # a function or callable object that takes a time and a number of
        # samples and returns values looked up at that time.
        get_samples,
        # the minimum pitch shift factor (simply used to limit the rate signal)
        ps_min=0.25,
        # the maxmimum pitch shift factor (used to limit the rate signal and
        # allocate a large enough ring buffer)
        ps_max=4.,
        # the block size of processing. This size is what is passed to get_samples.
        B=256,
        # Get a function to do the interpolation of the signal values
        get_interpolator=default_get_interpolator,
        # Given a minimum and maximum value that you want to look up, get the
        # range of values the interpolator needs in order to function.
        # must return integers
        get_interpolator_range=default_get_interpolator_range,
        # Given a number of unique values that need to be looked up, return the
        # number of values the interpolator needs in order to give all the
        # values. This works because the interpolator will be looking up values
        # in the range [0,N-1], and so naturally needs some extra points outside
        # of this range in order to do non-trivial interpolation. E.g., for
        # cubic Lagrange interpolation, you will [-1,N,N+1] or N+3 points.
        # This is used to comput the size of the ringbuffer.
        # must return an integer
        get_interpolator_n_points=default_get_interpolator_n_points,
        dtype=np.float64):

        # The most values that could be requested are
        # ceil((B*ps_max+get_interpolator_n_points(B))/B)*B
        tmp=B*ps_max+get_interpolator_n_points(B)
        sig_rb_size=0
        while sig_rb_size < tmp:
            sig_rb_size += B
        self.sig_rb = ringbuffer(sig_rb_size,dtype)

        # These values are invalid until after the first call to process
        # This is the index of the value at index 0 in the ring buffer
        self.sig_rb_idcs_valid = False
        self.sig_rb_min_idx = -1
        # This is the index of the value at index
        # self.sig_rb.contents_size() - 1 in the ringbuffer
        self.sig_rb_max_idx = -1

        self.get_samples = get_samples
        self.ps_min = ps_min
        self.ps_max = ps_max
        self.B = B
        self.get_interpolator=get_interpolator
        self.get_interpolator_range=get_interpolator_range
        self.dtype=dtype

        self.time_at_block_start=0
        self.pos_at_block_start=0

    def process(self,rate_sig):
        """
        rate_sig is an array of B values representing the time-increment. These
        values will be restricted between self.ps_min and self.ps_max. These
        values are accumulated to give the signal
        pos_signal=[self.pos_at_block_start,
        self.pos_at_block_start+rate_sig[0],
        ...,
        self.pos_at_block_start+sum(rate_sig)].
        The values pos_signal[:B] are used to look up values (using
        interpolation) and self.time_at_block_start is set to pos_signal[B].
        """
        # restrict rate signal
        rate_sig[rate_sig>self.ps_max] = self.ps_max
        rate_sig[rate_sig<self.ps_min] = self.ps_min

        # calculate the position signal by accumulating the rate signal
        pos_signal=np.zeros(self.B+1,dtype=self.dtype)
        np.cumsum(rate_sig,out=pos_signal[1:])
        pos_signal += self.pos_at_block_start
        first_required_idx,last_required_idx=self.get_interpolator_range(
            pos_signal[0],pos_signal[-2])

        if self.sig_rb_idcs_valid:
            fetch_start_idx = self.sig_rb_max_idx + 1
            # discard all values between self.sig_rb_min_idx and the first_required_idx
            n_discard = first_required_idx - self.sig_rb_min_idx
            self.sig_rb.advance_head(n_discard)
        else:
            fetch_start_idx = first_required_idx
            # start here so the accumlated value is correct (see the while loop below)
            self.sig_rb_max_idx = fetch_start_idx - 1

        self.sig_rb_min_idx = first_required_idx

        n_fetch_vals = last_required_idx - fetch_start_idx + 1

        # linear interpolation giving the look up times from the signal indices
        fetch_time_interpolator=interpolate.interp1d(
            [pos_signal[0],pos_signal[-1]],
            [self.time_at_block_start,self.time_at_block_start+self.B],
            fill_value='extrapolate')

        # fetch the required values
        while fetch_start_idx <= last_required_idx:
            fetch_time = fetch_time_interpolator(fetch_start_idx)
            x=self.get_samples(fetch_time,self.B)
            self.sig_rb.push_copy(x)
            fetch_start_idx += self.B
            self.sig_rb_max_idx += self.B
        self.sig_rb_idcs_valid = True
        
        # get the signal result by interpolating
        signal_interpolator=self.get_interpolator(
            np.arange(first_required_idx,last_required_idx+1,1),
            self.sig_rb.get_region(0,last_required_idx-first_required_idx+1))
        y=signal_interpolator(pos_signal[:-1])
        self.time_at_block_start += self.B
        self.pos_at_block_start = pos_signal[-1]
        return y
