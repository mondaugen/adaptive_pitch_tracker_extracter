/* Estimate the time of attacks in audio signals */

#include <math.h>
#include <string.h>
#include "pvoc_windows.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define SWAP(a,b)\
    ({ typeof(a) tmp = a;\
       a = b;\
       b = tmp; })

struct spec_diff_result {
    float *spec_diff;
    unsigned int length;
};

static void
spec_diff_result_free(struct spec_diff_result *x)
{
    if (x) {
        if (x->spec_diff) { free(x->spec_diff); }
        free(x);
    }
}

/* TODO: The policy should be that structs whose fields are known (i.e., users
of the struct can allocate them and access their fields without helper functions
provided by the implementation), should be freeable with free */
static spec_diff_result_new(unsigned int length)
{
    struct spec_diff_result *ret = calloc(sizeof(struct spec_diff_result),1);
    if (!ret) { goto fail; }
    ret->spec_diff = calloc(sizeof(float),length);
    if (!ret->spec_diff) { goto fail; }
    ret->length = length;
    return ret;
fail:
    spec_diff_result_free(ret);
    return NULL;
}

struct spec_diff_finder {
    unsigned int H;
    unsigned int W;
    const float *window;
    /* Auxiliary data for DFT calculations */
    void *dft_aux;
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
    if (init->init_window(ret->window,init->init_window_aux,init-W)) {
        goto fail;
    };
    ret->dft_aux = spec_diff_dft_new(init->W);
    if (!ret->dft_aux) { goto fail; }
    return ret;
fail:
    spec_diff_finder_free(ret);
    return NULL;
}

static struct spec_diff_result *
spec_diff_finder_find(struct spec_diff_finder *finder,
                      const float *x,
                      unsigned int len_x)
{
    /* The length of the array of complex values e.g., for an "unpacked"
    real-only DFT, this length is floor(finder->W/2) + 1 */
    unsigned int vz_length = spec_diff_complex_size(finder->W),
                 L = MAX(vz_length,finder->W);
    char tmp[L*sizeof(float complex)];
    float complex *z_tmp = (float complex *)tmp;
    float _X0[L], _X1[L],
          *X0 = _X0, *X1 = _X1,
          x_tmp = (float*)tmp;
    if (len_x == 0) { return NULL; }
    memset(X1,0,sizeof(float)*L);
    unsigned int length_spec_diff = pvoc_calc_n_blocks(
                                    len_x, finder->H, finder->W),
                 h,
                 n = 0;
    struct spec_diff_result *ret = spec_diff_result_new(length_spec_diff);
    if (!ret) { return NULL; }
    for (h = 0; h < length_spec_diff; h++) {
        /* multiply by window */
        dspm_mul_vf32_vf32_vf32(x + n,
                                finder->window,
                                X0,
                                finder->W);
        /* compute forward DFT of signal-window product */
        spec_diff_forward_dft(finder->dft_aux,X0,z_tmp);
        /* Compute absolute values of spectrum, real part will then contain the
        magnitude of the complex number and the imaginary part will be 0 */
        dspm_abs_vz32_vf32(z_tmp, X0, vz_length);
        /* compute difference between this frame and last */
        dspm_sub_vf32_vf32_vf32(X0,X1,x_tmp,vz_length);
        /* clip to 0 (we only want positive spectral flux) */
        dspm_clip_below_vf32_f32(x_tmp,0,vz_length);
        /* sum to get spectral difference */
        ret->spec_diff[h] = dspm_sum_vf32(x_tmp,vz_length);
        /* Make X0 into X1 and X1 is now X0, which will be overwritten on next
        pass */
        SWAP(X0,X1);
    }
    return ret;
}

struct attacks_from_spec_diff_finder {
    struct spec_diff_finder *spec_diff_finder;
    /* hop size */
    unsigned int H;
    /* window size */
    unsigned int W;
    /* smoothing factor */
    float smoothing;
    /* the local maxmimum filter coefficient computed from lmax_filt_rate*/
    float lmfr;
    /* the linear threshold in dB */
    float ng_th_A;
};

void
attacks_from_spec_diff_finder_free(struct attacks_from_spec_diff_finder *finder)
{
    if (finder) {
        if (finder->spec_diff_finder) {
            spec_diff_finder_free(finder->spec_diff_finder);
        }
        free(finder);
    }
}

static int get_normalized_window (float *w,
                                   void *_type,
                                   unsigned int length)
{
    char *type = _type;
    int ret = get_pvoc_window(w,type,length);
    if (ret) { return ret; }
    float win_sum = dspm_sum_vf32(w, length);
    /* TODO: Guarantee this error code be different from get_pvoc_window? */
    if (win_sum == 0) { return -2; }
    dspm_div_vf32_f32(w,win_sum,length);
    return 0;
}

static float 
compute_lmfr(float lmax_filt_rate, float H)
{
    /*
    calculate the max filter rate
    number of hops in the lmax_filt_rate
    */
    float lmfr_n_H=args->lmax_filt_rate/H;
    /* Must be greater than 0 */
    if (lmfr_n_H <= 0) { lmfr_n_H=1; }
    /* what number to this power is .01 ? */
    float ret = pow(.01,1/lmfr_n_H);
    return ret;
}

static int
attacks_from_spec_diff_finder_chk_args(
struct attacks_from_spec_diff_finder_args *args)
{
    if (args->H < 1) { return -1; }
    if (args->W < 1) { return -2; }
    if (!args->window_type) { return -3; }
    if ((args->smoothing < 0) || (args->smoothing > 1)) { return -4; }
    if (args->lmax_filt_rate < 1) { return -5; }
    return 0;
}

struct attacks_from_spec_diff_finder *
attacks_from_spec_diff_finder_new(
struct attacks_from_spec_diff_finder_args *args)
{
    if(attacks_from_spec_diff_finder_chk_args(args)) { return NULL; }
    struct attacks_from_spec_diff_finder *ret = calloc(1,
        sizeof(struct attacks_from_spec_diff_finder));
    if (!ret) { return NULL; }
    struct spec_diff_finder_init spec_diff_finder_init = {
        .H = args->H,
        .W = args->W
        .init_window = get_normalized_window,
        .init_window_aux = args->window_type
    };
    ret->spec_diff_finder = spec_diff_finder_new(&spec_diff_finder_init);
    if (!ret->spec_diff_finder) { goto fail; }
    ret->H = args->H;
    ret->W = args->W;
    ret->smoothing = smoothing;
    ret->lmfr = compute_lmfr(args->lmax_filt_rate, args->H);
    return ret;
fail:
    attacks_from_spec_diff_finder_free(ret);
    return NULL;
}

/*
smooth values using an IIR filter
# a = 1 uses 100% of current value in output and 0% of past values
# a = 0.5 uses 50% current value, and 50% past values, etc.
y=signal.lfilter([a],[1,-(1-a)],x)
*/
static void
iir_avg(const float *x, float *y, unsigned int length, float a)
{
    float y_1 = 0, y0;
    while (length--) {
        y0 = a * *x++ + (1 - a) * y_1;
        *y++ = y0;
        y_1 = y0;
    }
}

/* Returns pairs that are an estimation of the attack start and end times. */
struct attacks_from_spec_diff_result *
attacks_from_spec_diff_finder_compute(
    struct attacks_from_spec_diff_finder *finder,
    float *x
    unsigned int length)
{
    if (!x) { return NULL; }
    if (length < 1) { return NULL; }
    struct spec_diff_result *sd = NULL;
    struct index_points *sd_mins = NULL,
                        *sd_maxs = NULL,
                        *sd_mins_filtered = NULL;
    struct attacks_from_spec_diff_result *ret = NULL;
    /* Find spectral difference */
    sd = spec_diff_finder_find(finder->spec_diff_finder,x,length);
    if (!sd) { goto fail; }
    /* Smooth the spectral difference */
    iir_avg(sd->spec_diff,sd->spec_diff,sd->length,finder->smoothing);
fail:
    if (sd) { spec_diff_result_free(sd); }
    if (sd_mins) { index_points_free(sd_mins); }
    if (sd_maxs) { index_points_free(sd_maxs); }
    if (sd_mins_filtered) { index_points_free(sd_mins_filtered); }
    return ret;
}

                                      
