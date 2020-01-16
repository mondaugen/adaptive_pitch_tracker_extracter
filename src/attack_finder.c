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

struct spec_diff_finder {
    unsigned int H;
    unsigned int W;
    const float *window;
    /* Auxiliary data for DFT calculations */
    void *dft_aux;
    /* This stores the complex data and eventually the magnitude data. This can
    be anything but the first sizeof(float)*W bytes must be floats and at some
    point will represent the magnitude spectrum, which will be used to calculate
    the spectral difference. */
    void *X0;
    /* This is the same as X0 but stores the data from a hop ago. This needs to
    be initialized to contain 0s for the first sizeof(float)*W bytes. */
    void *X_1;
};

struct spec_diff_finder_init {
    unsigned int H;
    unsigned int W;
    void (*init_window)(float *w, void *aux, unsigned int window_length);
    void *init_window_aux;
};

static int
spec_diff_finder_init_chk_args(struct spec_diff_finder_init *init)
{
    if (init->H == 0) { return -1; }
    if (init->W == 0) { return -2; }
    if (!init->init_window) { return -3; }
    return 0;
}

void
spec_diff_finder_free(struct spec_diff_finder *finder)
{
    if (finder) {
        if (finder->window) { free(finder->window); }
        if (finder->dft_aux) { spec_diff_dft_free(finder->dft_aux); }
        if (finder->X0) { spec_diff_complex_free(finder->X0); }
        if (finder->X1) { spec_diff_complex_free(finder->X1); }
        free(finder);
    }
}

struct spec_diff_finder *
spec_diff_finder_new(struct spec_diff_finder_init *init)
{
    if (spec_diff_finder_init_chk_args(init)) { return NULL; }
    struct spec_diff_finder *ret = calloc(1,sizeof(struct spec_diff_finder));
    if (!ret) { goto fail; }
    ret->H = init->H;
    ret->W = init->W;
    ret->window = calloc(init->W,sizeof(float));
    if (!ret->window) { goto fail; }
    ret->dft_aux = spec_diff_dft_new(init->W);
    if (!ret->dft_aux) { goto fail; }
    ret->X0 = spec_diff_complex_alloc(init->W);
    if (!ret->X0) { goto fail; }
    ret->X1 = spec_diff_complex_alloc(init->W);
    if (!ret->X1) { goto fail; }
    return ret;
fail:
    spec_diff_finder_free(ret);
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
