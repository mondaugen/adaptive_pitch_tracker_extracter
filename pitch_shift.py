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
        # a function or callable object that takes a time in samples and returns
        # B samples looked up at that time.
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
        # This is used to compute the size of the ringbuffer and must return an
        # integer.
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

    def set_pos_at_block_start(self,pos):
        self.pos_at_block_start=pos
        self.time_at_block_start=pos
        self.sig_rb_idcs_valid=False

    def process_pos_sig(self,ps_pos_sig,ts_pos_sig):
        # TODO In the C implementation, ps_pos_sig should have at least 16 bits
        # of fractional precision to keep the tuning error less than 1% of a
        # cent

        B=len(ps_pos_sig)-1

        first_required_idx,last_required_idx=self.get_interpolator_range(
            ps_pos_sig[0],ps_pos_sig[-2])

        if self.sig_rb_idcs_valid:
            fetch_start_idx = self.sig_rb_max_idx + 1
            # discard all values between self.sig_rb_min_idx and the first_required_idx
            n_discard = first_required_idx - self.sig_rb_min_idx
            self.sig_rb.advance_head(n_discard)
        else:
            self.sig_rb.reset()
            fetch_start_idx = first_required_idx
            # start here so the accumlated value is correct (see the while loop below)
            self.sig_rb_max_idx = fetch_start_idx - 1

        self.sig_rb_min_idx = first_required_idx

        n_fetch_vals = last_required_idx - fetch_start_idx + 1

        # linear interpolation giving the look up times from the signal indices
        fetch_time_interpolator=interpolate.interp1d(
            [ps_pos_sig[0],ps_pos_sig[-1]],
            [ts_pos_sig[0],ts_pos_sig[-1]],
            fill_value='extrapolate')

        # fetch the required values
        while fetch_start_idx <= last_required_idx:
            fetch_time = fetch_time_interpolator(fetch_start_idx)
            x=self.get_samples(fetch_time)
            self.sig_rb.push_copy(x)
            fetch_start_idx += self.B
            self.sig_rb_max_idx += self.B
        # TODO why is this made true here?
        self.sig_rb_idcs_valid = True
        
        # get the signal result by interpolating
        signal_interpolator=self.get_interpolator(
            np.arange(first_required_idx,last_required_idx+1,1),
            self.sig_rb.get_region(0,last_required_idx-first_required_idx+1))
        y=signal_interpolator(ps_pos_sig[:-1])
        self.time_at_block_start = ts_pos_sig[-1]
        self.pos_at_block_start = ps_pos_sig[-1]
        return y

    def process(self,ps_rate_sig,ts_rate_sig=None):
        """
        ps_rate_sig is an array of B values representing the time-increment for
        time-stretching. These values will be restricted between self.ps_min and
        self.ps_max. These values are accumulated to give the signal
        ps_pos_sig=[self.pos_at_block_start,
        self.pos_at_block_start+ps_rate_sig[0],
        ...,
        self.pos_at_block_start+sum(ps_rate_sig)].
        The values ps_pos_sig[:B] are used to look up values (using
        interpolation) and self.time_at_block_start is set to ps_pos_sig[B].
        It is confusing, but here the ps_rate_sig and the accumulated ps_pos_sig
        represent the pitch shift amount, which will be counteracted by time
        stretching at an inverse rate.
        """

        if ts_rate_sig is None:
            ts_rate_sig=np.zeros_like(ps_rate_sig)
            ts_rate_sig[:]=1

        # B is length of ps_rate_sig
        B=len(ps_rate_sig)

        # restrict rate signal
        ps_rate_sig[ps_rate_sig>self.ps_max] = self.ps_max
        ps_rate_sig[ps_rate_sig<self.ps_min] = self.ps_min

        # calculate the position signal by accumulating the rate signal
        # This ps_pos_sig is the one formed by summing the ps_rate_sig only
        # (not after multiplying by the time_stretch signal)
        ps_pos_sig=np.zeros(B+1,dtype=self.dtype)
        np.cumsum(ps_rate_sig,out=ps_pos_sig[1:])

        ts_pos_sig=np.zeros(B+1,dtype=self.dtype)
        np.cumsum(ts_rate_sig,out=ts_pos_sig[1:])

        # TODO With the u24q8 format, we can address an array of length of 2**24
        # - get_interpolator_n_points(0) (about 350 seconds at 48KHz sample
        # rate).  This length is reasonable for our purposes: we assume that the
        # get_samples implementation reduces the requested index modulo the
        # length of the sound-array. However also after about 350 seconds
        # (depending on the pitch-shift factor) the ps_pos_sig will roll-over.
        # It is not unreasonable that this pitch-shifter might be playing
        # through a sound (looping) for more than 5 minutes, so we need to also
        # reduce these positions modulo the sound-array length.  In the C
        # implementation, we will have ps_pos_sig be in u48q16 format and
        # ts_pos_sig in s48q16 format. Before interpolation takes place we
        # subtract from the data-point x-positions and the looked up x positions
        # to put the values in a reasonable range whose maximum index fits in
        # u24q8. This slight increase in computation makes it so don't need to
        # worry about roll over for almost a century
        ps_pos_sig += self.pos_at_block_start
        ts_pos_sig += self.time_at_block_start

        return self.process_pos_sig(ps_pos_sig,ts_pos_sig)

