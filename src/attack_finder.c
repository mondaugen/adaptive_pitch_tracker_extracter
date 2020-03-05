/* Estimate the time of attacks in audio signals */

#include <stdint.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "pvoc_windows.h"
#include "attack_finder.h"
#include "dsp_math.h"
#include "find_extrema.h"

#define MAX(a,b) \
  ({ typeof (a) _a = (a); \
     typeof (b) _b = (b); \
     _a > _b ? _a : _b; })

#define SWAP(a,b)\
    ({ typeof(a) tmp = a;\
       a = b;\
       b = tmp; })

/*
Define these somewhere for the relevant implementation.
The argument for doing it this way is because there are many ways to
implement the DFT. For example you might have a real-only DFT that is packed or
unpacked, or one where the values are permuted. These may all exist one day in
the dsp_math routines, but only one can be chosen for this routine.
*/
/* free a dft object */
extern void spec_diff_dft_free(void *aux);
/* allocate a dft object for a given window size W */
extern void *spec_diff_dft_new(unsigned int W);
/* carry out a forward DFT */
extern void spec_diff_forward_dft(void *aux, float *time, float complex *freq);
/* Return the number of complex values in the result of a DFT from W real
numbers. For example if the DFT is the unpacked real DFT then the number of
complex values is floor(W/2) + 1 for most implementations. */
extern unsigned int spec_diff_complex_size(unsigned int W);

struct spec_diff_result {
    float *spec_diff;
    unsigned int length;
};

struct spec_diff_finder_init {
    unsigned int H;
    unsigned int W;
    /* Returns non-zero if it didn't succeed in initializing the window. NOTE:
    From prototype implementation in spectral_difference.py, the window W is
    divided by sum(W) to normalize it, so the init_window should also do this
    normalization to the window it puts in w. */
    int (*init_window)(float *w, void *aux, unsigned int window_length);
    void *init_window_aux;
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
static struct spec_diff_result *
spec_diff_result_new(unsigned int length)
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

unsigned int *
attacks_result_extract_beginnings(
    struct attacks_from_spec_diff_result *a)
{
    unsigned int *ret = malloc(sizeof(unsigned int)*a->n_attack_time_pairs),
                  n;
    if (!ret) { return NULL; }
    for (n = 0; n < a->n_attack_time_pairs; n++) {
        ret[n] = a->attack_time_pairs[n].beginning;
    }
    return ret;
}

unsigned int *
attacks_result_extract_ends(
    struct attacks_from_spec_diff_result *a)
{
    unsigned int *ret = malloc(sizeof(unsigned int)*a->n_attack_time_pairs),
                  n;
    if (!ret) { return NULL; }
    for (n = 0; n < a->n_attack_time_pairs; n++) {
        ret[n] = a->attack_time_pairs[n].end;
    }
    return ret;
}

struct spec_diff_finder {
    unsigned int H;
    unsigned int W;
    float *window;
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
    if (init->init_window(ret->window,init->init_window_aux,init->W)) {
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
          *x_tmp = (float*)tmp;
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
        n += finder->H;
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
    /* the linear threshold */
    float ng_th_A;
    /* the spectral difference threshold */
    float sd_th;
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
    float lmfr_n_H=lmax_filt_rate/H;
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
        .W = args->W,
        .init_window = get_normalized_window,
        /* We can discard const because get_normalized_window only reads init_window_aux */
        .init_window_aux = (void*)(uintptr_t)args->window_type
    };
    ret->spec_diff_finder = spec_diff_finder_new(&spec_diff_finder_init);
    if (!ret->spec_diff_finder) { goto fail; }
    ret->H = args->H;
    ret->W = args->W;
    ret->smoothing = args->smoothing;
    ret->lmfr = compute_lmfr(args->lmax_filt_rate, args->H);
    /* To avoid having to find the mean square when computing, we scale the
    threshold accordingly */
    ret->ng_th_A = ret->W*pow(10,args->ng_th/10);
    ret->sd_th = args->sd_th;
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

static inline float *
local_sum_square(
    const float *x,
    unsigned int N,
    unsigned int H,
    unsigned int W)
{
    unsigned int n_mean_square = pvoc_calc_n_blocks(N, H, W),
                 n = 0;
    float tmp[W], *ret = NULL, *pret;
    const float *px = x;
    ret = malloc(n_mean_square*sizeof(float));
    if (!ret) { return NULL; }
    pret = ret;
    for (n = 0; n < n_mean_square; n++) {
        dspm_mul_vf32_vf32_vf32(px, px, tmp, W);
        *pret++ = dspm_sum_vf32(tmp,W);
        px += H;
    }
    return ret;
}

static inline float
compute_gate(
    const float *px,
    unsigned int W,
    float thresh)
{
    float tmp[W], rms;
    dspm_mul_vf32_vf32_vf32(px,px,tmp,W);
    rms = dspm_sum_vf32(tmp,W);
    if (rms > thresh) { return 1.; }
    return 0;
}

/*
What is returned is the indices of the gate changes. The gate is either 0 or
1 and has the initial value 0. The first index is when the gate goes from 0 to
1, then the second when it goes from 1 to 0, etc.
*/
static inline unsigned int *
gate_via_rms_thresh(
    /* signal to determine gates from */
    const float *x,
    /* its length */
    unsigned int len_x,
    /* hop size */
    unsigned int H,
    /* window size */
    unsigned int W,
    /* after returns contains the length of the returned array */
    unsigned int *n_gate_changes,
    /* if RMS of x is greater than this threshold, gate is high */
    float thresh)
{
    /* algorithm done in 2 passes to save on memory */
    unsigned int *gates = NULL,
                 n_mean_square = pvoc_calc_n_blocks(len_x, H, W),
                 n, ngc, pass;
    float gate, gate_1, dgate;
    const float *px;
    for (pass = 0; pass < 2; pass++) {
        px = x;
        gate_1 = 0;
        ngc = 0;
        for (n = 0; n < n_mean_square; n++) {
            gate = compute_gate(px,W,thresh);
            dgate = gate - gate_1;
            if (dgate) { 
                /* first pass we just count the number of gates */
                if (pass > 0) { gates[ngc] = n; }
                ngc++;
            }
            gate_1 = gate;
            px += H;
        }
        if (pass > 0) { break; }
        gates = malloc(sizeof(unsigned int)*ngc);
        if (!gates) { return NULL; }
    }
    *n_gate_changes = ngc;
    return gates;
}

/*
filter out the indices in find_closest by leaving only the ones coming right
after indices on the left.
Output is sorted in ascending order.
filter and find_closest must be sorted in ascending order and contain only unique items
*/
unsigned int *
attack_finder_closest_index_after(
    const unsigned int *filtered,
    /* length of filtered */
    unsigned int n_filtered,
    const unsigned int *find_closest,
    /* length of find_closest */
    unsigned int n_find_closest,
    /* after returning, contains the number of indices in resulting array */
    unsigned int *n_idcs,
    /* if non-zero, filters out indices in find_closest by leaving only those
    coming before indices on the right */ 
    int reverse)
{
    /* TODO: Is this the best policy? This function shouldn't get called at all
    in the following cases */
    if ((!filtered)||
        (n_filtered==0)||
        (!find_closest)||
        (n_find_closest==0)||
        (!n_idcs)) {
        return NULL;
    }
    unsigned int ary_len = MAX(filtered[n_filtered-1], 
                 find_closest[n_find_closest-1]) + 1,
                 *idcs = NULL,
                 n, n_idcs_, val, pass, idc, f_idc, fc_idc;
    float x_1, x, dx, a_n, b_n;
    for (pass = 0; pass < 2; pass++) {
        x_1 = 0;
        n_idcs_ = 0;
        f_idc = 0;
        fc_idc = 0;
        for (n = 0; n < ary_len; n++) {
            a_n = 0; b_n = 0;
            if (f_idc < n_filtered) {
                idc = reverse ? n_filtered - f_idc - 1 : f_idc;
                val = filtered[idc];
                if (reverse) { val = ary_len - val - 1; }
                if (val == n) {
                    a_n = 1;
                    f_idc++;
                }
            }
            if (fc_idc < n_find_closest) {
                idc = reverse ? n_find_closest - fc_idc - 1 : fc_idc;
                val = find_closest[idc];
                if (reverse) { val = ary_len - val - 1; }
                if (val == n) {
                    b_n = 1;
                    fc_idc++;
                }
            }
            x = x_1 + a_n - b_n;
            if (x < 0) { x = 0; }
            if (x > 1) { x = 1; }
            dx = x - x_1;
            x_1 = x;
            if (dx < 0) {
                val = n;
                if (reverse) { val = ary_len - val - 1; }
                if (pass > 0) {
                    idcs[n_idcs_] = val;
                }
                n_idcs_++;
            }
        }
        if (pass == 0) {
            idcs = malloc(n_idcs_*sizeof(unsigned int));
            if (!idcs) { return NULL; }
        }
    }
    *n_idcs = n_idcs_;
    if (reverse) {
        /* the indices are currently in decreasing order, reverse the order so
        they are ascending */
        dspm_rev_vu32(idcs,n_idcs_);
    }
    return idcs;
}

static inline void threshold_to_gate(
    float *x,
    unsigned int N,
    float threshold)
{
    while (N--) {
        *x = *x > threshold ? 1 : 0;
        x++;
    }
}

/* 
idcs is an array of length *n_idcs containing unique unsigned ints in ascending
order.
gate_changes is an array of length n_gate_changes containing unique unsigned
integers in ascending order that indicate when the gate goes from 0 to 1 or 1 to
0.
len_sig is the length of the signal gate_changes is derived from, i.e., neither
idcs nor gate_changes should contain an index >= len_sig
*/
unsigned int *
attack_finder_index_mask(
    /* gets freed by this function if it succeeds */
    unsigned int *idcs,
    /* on entry contains length of indcs,
    on exit contains length of masked idcs */
    unsigned int *n_idcs,
    /* all values in gate_changes must be unique 
    and sorted in ascending order */
    const unsigned int *gate_changes,
    /* assumes n_gate >= max(idcs) */
    unsigned int n_gate_changes,
    /* length of signal that gate is based off of */
    unsigned int len_sig)
{
    unsigned int n, n_entry_idcs = *n_idcs,
                 n_filtered, *ret = NULL, pass, gcn, n_idc;
    int gate, gate_dir;
    for (pass = 0; pass < 2; pass++) {
        gcn = 0;
        gate = 0;
        gate_dir = 1;
        n_idc = 0;
        n_filtered = 0;
        for (n = 0; n < len_sig; n++) {
            if ((gcn < n_gate_changes) &&
                    (gate_changes[gcn] == n)) {
                gate += gate_dir;
                gate_dir *= -1;
                gcn++;
            }
            if ((idcs[n_idc] == n)) {
                if (gate) {
                    if (pass > 0) { ret[n_filtered] = idcs[n_idc]; }
                    n_filtered++;
                    /* force quit if no more indices */
                    if (n_idc >= n_entry_idcs) { n = len_sig; }
                }
                n_idc++;
            }
        }
        if (pass > 0) { break; }
        ret = malloc(sizeof(unsigned int)*n_filtered);
        /* If failed just return idcs */
        if (!ret) { ret = idcs; goto fail; }
    }
    free(idcs);
    *n_idcs = n_filtered;
fail:
    return ret;
}

void
attacks_from_spec_diff_result_free(struct attacks_from_spec_diff_result *r)
{
    if (r) {
        if (r->attack_time_pairs) {
            free(r->attack_time_pairs);
        }
        free(r);
    }
}

static struct attacks_from_spec_diff_result *
attacks_from_spec_diff_result_new(unsigned int n_attack_time_pairs)
{
    struct attacks_from_spec_diff_result *ret = NULL;
    ret = calloc(1,sizeof(struct attacks_from_spec_diff_result));
    if (!ret) { goto fail; }
    ret->attack_time_pairs = calloc(n_attack_time_pairs,
        sizeof(struct attack_sample_time_pair));
    if (!ret->attack_time_pairs) { goto fail; }
    return ret;
fail:
    attacks_from_spec_diff_result_free(ret);
    return NULL;
}

static struct attacks_from_spec_diff_result *
unpaired_ends(
    unsigned int *down,
    unsigned int n_down)
{
    struct attacks_from_spec_diff_result *ret =
        attacks_from_spec_diff_result_new(n_down);
    if (!ret) { return NULL; }
    ret->n_attack_time_pairs = 0;
    while (n_down--) {
        ret->attack_time_pairs[ret->n_attack_time_pairs++].end =
            *down++;
    }
    return ret;
}

/* up, down cannot be NULL */
static struct attacks_from_spec_diff_result *
up_down_match(
    unsigned int *up,
    unsigned int n_up,
    unsigned int *down,
    unsigned int n_down)
{
    struct attacks_from_spec_diff_result *ret = NULL;
    unsigned int n, m, pass, n_pairs;
    for (pass = 0; pass < 2; pass++) {
        n_pairs = 0;
        m = 0;
        for (n = 0; n < n_up; n++) {
            while ((m < n_down) && (down[m] < up[n])) {
                m++;
            }
            if (m < n_down) {
                if (pass > 0) {
                    ret->attack_time_pairs[n_pairs].beginning = up[n];
                    ret->attack_time_pairs[n_pairs].end = down[m];
                }
                n_pairs++;
            }
        }
        if (pass == 0) {
            ret = attacks_from_spec_diff_result_new(n_pairs);
            if (!ret) { return NULL; }
        }
    }
    ret->n_attack_time_pairs = n_pairs;
    return ret;
}
    
static void
offset_extrema(
    unsigned int *idcs,
    unsigned int n_idcs,
    struct attacks_from_spec_diff_finder *finder)
{
    /* add one hop to compensate for the differencing operation */
    dspm_add_vu32_u32(idcs,n_idcs,1);
    /* multiply by hop size because we want the indices in samples not hops */
    dspm_mul_vu32_u32(idcs,n_idcs,finder->H);
    /* compensate for window by offseting the maxima by half a window */
    dspm_add_vu32_u32(idcs,n_idcs,finder->W/2);
}

/* Returns pairs that are an estimation of the attack start and end times. */
struct attacks_from_spec_diff_result *
attacks_from_spec_diff_finder_compute(
    struct attacks_from_spec_diff_finder *finder,
    float *x,
    unsigned int length,
    /* If non-zero, trys to find the beginning of the attack just before the
    estimated attack peak and if found, puts it in the "beginning" field of the
    same attack_sample_time_pair as the attack peak time. Otherwise it just puts
    the attack time peaks in the "end" field of the attack_sample_time_pairs. */
    int return_time_pairs)
{
    if (!x) { return NULL; }
    if (length < 1) { return NULL; }
    struct spec_diff_result *sd = NULL;
    unsigned int *sd_mins = NULL,
                 *sd_maxs = NULL,
                 *sd_mins_filtered = NULL,
                 *sd_gate = NULL,
                 n_sd_maxs,
                 n_sd_mins,
                 n_sd_mins_filtered,
                 n_gate_changes;
    struct attacks_from_spec_diff_result *ret = NULL;
    float *thresh = NULL;
    /* Find spectral difference */
    sd = spec_diff_finder_find(finder->spec_diff_finder,x,length);
    if (!sd) { goto fail; }
    /* Smooth the spectral difference */
    iir_avg(sd->spec_diff,sd->spec_diff,sd->length,finder->smoothing);
    /* Find local maxima using a discounting technique */
#ifdef DEBUG
    thresh = calloc(sizeof(float),sd->length);
    if (!thresh) { goto fail; }
#endif
    sd_maxs = discount_local_max_f32(sd->spec_diff, sd->length, &n_sd_maxs,
        local_max_type_right, finder->lmfr, finder->sd_th,thresh);
    if (!sd_maxs || (n_sd_maxs == 0)) { goto fail; }
    sd_gate = gate_via_rms_thresh(x,length,finder->H,finder->W,&n_gate_changes,
    finder->ng_th_A);
    if (!sd_gate || (n_gate_changes == 0)) { goto fail; }
    sd_maxs = attack_finder_index_mask(
        sd_maxs,
        &n_sd_maxs,
        sd_gate,
        n_gate_changes,
        sd->length);
    if (!sd_maxs || (n_sd_maxs == 0)) { goto fail; }
    if (return_time_pairs) {
        /* Multiply spectral difference by -1 to find local minima */
        dspm_neg_vf32(sd->spec_diff, sd->length);
        /* Find local minima */
        sd_mins = local_max_f32(sd->spec_diff, sd->length, &n_sd_mins,
            local_max_type_left);
        if (!sd_mins || (n_sd_mins == 0)) { goto fail; }
        sd_mins_filtered = attack_finder_closest_index_after(
            sd_maxs,
            n_sd_maxs,
            sd_mins,
            n_sd_mins,
            &n_sd_mins_filtered,
            1 /* reverse */);
        if (!sd_mins_filtered || (n_sd_mins_filtered == 0)) { goto fail; }
        offset_extrema(sd_maxs,n_sd_maxs,finder);
        offset_extrema(sd_mins_filtered,n_sd_mins_filtered,finder);
        ret = up_down_match(
            sd_mins_filtered,
            n_sd_mins_filtered,
            sd_maxs,
            n_sd_maxs);
    } else {
        offset_extrema(sd_maxs,n_sd_maxs,finder);
        ret = unpaired_ends(sd_maxs,n_sd_maxs);
    }
fail:
    if (sd) { spec_diff_result_free(sd); }
    if (sd_mins) { free(sd_mins); }
    if (sd_maxs) { free(sd_maxs); }
    if (sd_mins_filtered) { free(sd_mins_filtered); }
    if (sd_gate) { free(sd_gate); }
    if (thresh) { free(thresh); }
    return ret;
}
