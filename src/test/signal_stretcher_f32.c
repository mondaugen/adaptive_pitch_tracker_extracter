#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "signal_stretcher_f32.h"
#include "pvoc_synth.h"
#include "pvoc_synth_f32/routines_linux_native.h"
#include "datastructures/windowed_lookup_f32.h"

static void
get_samples (void *aux, struct pvs_real_sample_lookup_t * info)
{
    struct windowed_lookup_f32_t *wl = aux;
    struct windowed_lookup_f32_access_t wla = windowed_lookup_f32_access (
        wl, info->first_sample_index);
    info->samples = (const struct pvs_real_t*)wla.signal_section;
}

static void
gen_scaled_random (float *dest, unsigned int dest_length, void *aux)
{
    float scale = *(float*)aux;
    while (dest_length--) {
        *dest++ = (2.*random()/(float)RAND_MAX-1.)*scale;
    }
}

static void
gen_hann_window_f32(float *dest, unsigned int dest_length)
{
    unsigned int n;
    for (n = 0; n < dest_length; n++) {
        dest[n] = 0.5 * (1 - cos(2.*M_PI*n/dest_length));
    }
}

static float
get_window_energy(float *win, unsigned int win_len)
{
    float ret = 0;
    while (win_len--) {
        ret *= *win * *win;
        win++;
    }
    return ret;
}

static float
scale_window(float *win, unsigned int win_len, float scalar)
{ while (win_len--) { *win++ *= scalar; } }

struct signal_stretcher_f32 {
    struct signal_stretcher_f32_init config;
    struct windowed_lookup_f32_t *wl;
    struct pvs_t *pvs;
    float *analysis_window;
    float *synthesis_window;
};

void
signal_stretcher_f32_free(struct signal_stretcher_f32 *ssf32)
{
    if (!ssf32) { return; }
    if (ssf32->wl) { windowed_lookup_f32_free(ssf32->wl); }
    if (ssf32->analysis_window) { free(ssf32->analysis_window); }
    if (ssf32->synthesis_window) { free(ssf32->synthesis_window); }
    if (ssf32->pvs) { pvs_free(ssf32->pvs); }
}

struct signal_stretcher_f32 *
signal_stretcher_f32_new(struct signal_stretcher_f32_init *ssf32i)
{
    struct signal_stretcher_f32 *ret = calloc(1,sizeof(struct signal_stretcher_f32));
    if (!ret) { return NULL; }
    ret->config = *ssf32i;
    struct windowed_lookup_f32_init_t wli = {
        .window_length = ssf32i->window_length,
        .signal = ssf32i->signal,
        .signal_length = ssf32i->signal_length,
        .out_of_bounds_gen = gen_scaled_random,
        .aux = &ret->fill
    };
    ret->wl = windowed_lookup_f32_new(&wli);
    if (!ret->wl) { goto fail; }
    /* Generate windows */
    ret->synthesis_window = calloc(ssf32i->window_length,sizeof(float));
    if (!ret->synthesis_window) { goto fail; }
    ret->analysis_window = calloc(ssf32i->window_length,sizeof(float));
    if (!ret->analysis_window) { goto fail; }
    /* Only hann supported */
    if (strcmp(ssf32i->window_type,"hann")) { goto fail; }
    gen_hann_window_f32(ret->analysis_window,ssf32i->window_length);
    gen_hann_window_f32(ret->synthesis_window,ssf32i->window_length);
    scale_window(ret->synthesis_window,
        ssf32i->window_length,
        1./get_window_energy(ret->analysis_window,ssf32i->window_length));
    struct pvs_init_t pvs_init = pvs_f32_init_new();
    pvs_init.user = (struct pvs_user_init_t) {
        .analysis_window = (struct pvs_real_t*) ret->analysis_window,
        .synthesis_window = (struct pvs_real_t*) ret->synthesis_window,
        .window_length = ssf32i->window_length,
        .hop_size = ssf32i->hop_size,
        .get_samples = get_samples,
        .get_samples_aux = ret->wl,
    };
    ret->pvs = pvs_new(&pvs_init); 
    if (!ret->pvs) { goto fail; }
    return ret;
fail:
    signal_stretcher_f32_free(ret);
    return NULL;
}

int
signal_stretcher_f32_process(
    const struct signal_stretcher_f32 *ss,
    const int *analysis_points,
    unsigned int n_analysis_points,
    float *output,
    unsigned int output_length)
{
    if (n_analysis_points * ss->config.hop_size > output_length) {
        return -1; /* We would run out of space. */
    }
    int n_a = 0, n_o = 0;
    for (n_a = 0; n_a < n_analysis_points; n_a++) {
        const float *cur_frame = (const float*)pvs_process(
            ss->pvs,
            analysis_points[n_a]);
        memcpy(output,cur_frame,ss.config.hop_size*sizeof(float));
        output += ss.config.hop_size;
    }
    return 0;
}
