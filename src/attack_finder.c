/* Estimate the time of attacks in audio signals */

#include "pvoc_windows.h"

struct spec_diff_result {
    float *spec_diff;
    unsigned int length;
};

static void
spec_diff_result_free(struct spec_diff_result *x)
{
    if (x->spec_diff) { free(x->spec_diff); }
    free(x);
}

static struct spec_diff_result *
spec_diff(const float *x,
          unsigned int len_x,
          unsigned int H,
          unsigned int W,
          const char *window_type)
{
    float w[W];
    /* TODO: Figure out how to determine spec_diff length */
    unsigned int length_spec_diff;
    if (get_pvoc_window(w,window_type,W)) {
        return NULL;
    }
/*
def spectral_diff(x,H,W,window_type):
    w=signal.get_window(window_type,W)
    x_framed=common.frame(x,H,W)*w[:,None]
    X_framed=np.fft.rfft(x_framed,axis=0)/np.sum(W)
    xd_abs=np.abs(X_framed)
    sd=xd_abs[:,1:]-xd_abs[:,:-1]
    sd[sd<0]=0
    sd=np.sum(sd,axis=0)
    return sd
*/

struct attacks_from_spec_diff_args {
    /* hop size */
    unsigned int H;
    /* window size */
    unsigned int W;
    /* window type */
    const char *window_type;
    /*
    smoothing factor s. 0 < s <= 1 1 means no smoothing, 0 would mean totally
    rejecting new values
    */
    float smoothing;
    /*
    max filter discount rate, the time in samples until it reaches 1% of the
    maximum
    */
    unsigned int lmax_filt_rate;
    /* the threshold in dB for the noise gate */
    float ng_th;
};

#define attacks_from_spec_diff_args_default \
(struct attacks_from_spec_diff_args) {\
    .H = 256,\
    .W = 1024,\
    .window_type = "hann",\
    .smoothing = 1,\
    .lmax_filt_rate = 16000,\
    ng_th=-60\
}

/* Used to store time of beginning and end of attack */
struct attack_sample_time_pair {
    unsigned int beginning;
    unsigned int end;
};

struct attacks_from_spec_diff_result {
    struct attack_sample_time_pair *attack_time_pairs;
    unsigned int n_attack_time_pairs;
};
