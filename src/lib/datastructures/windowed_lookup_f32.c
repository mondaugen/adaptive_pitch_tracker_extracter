/*
Data structure for lookups of f32 that will have a fixed window size and that
may request out of bounds data.
*/

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "windowed_lookup_f32.h"

struct windowed_lookup_f32_t {
    struct windowed_lookup_f32_init_t config;
    /* An array of length 2 * signal_length holding signal_length fill values
    followed by the first signal_length values */
    float *signal_start;
    /* An array of length 2 * signal_length holding the last signal_length
    values followed by signal_length fill values */
    float *signal_end;
};

static int init_arg_chk(struct windowed_lookup_f32_init_t *wli)
{
    if (!wli->signal) { return -1; }
    if (wli->window_length < 0) { return -2; }
    if (wli->signal_length < 0) { return -3; }
    if (wli->window_length > INT_MAX) { return -4; }
    if (wli->signal_length > INT_MAX) { return -5; }
    return 0;
}

void
windowed_lookup_f32_free(struct windowed_lookup_f32_t *wl)
{
    if (!wl) { return; }
    if (wl->signal_start) { free(wl->signal_start); }
    if (wl->signal_end) { free(wl->signal_end); }
    free(wl);
}

struct windowed_lookup_f32_t *
windowed_lookup_f32_new(struct windowed_lookup_f32_init_t *config)
{
    if (init_arg_chk(config)) { return NULL; }
    struct windowed_lookup_f32_t *ret = calloc(1,sizeof(struct windowed_lookup_f32_t));
    if (!ret) { goto fail; }
    ret->signal_start = calloc(2*config->window_length,sizeof(float));
    if (!ret->signal_start) { goto fail; }
    ret->signal_end = calloc(2*config->window_length,sizeof(float));
    if (!ret->signal_end) { goto fail; }
    if (config->out_of_bounds_gen) {
        config->out_of_bounds_gen(
            ret->signal_start,
            config->window_length,
            config->aux);
        config->out_of_bounds_gen(
            ret->signal_end+config->window_length,
            config->window_length,
            config->aux);
    }
    memcpy(
        ret->signal_start+config->window_length,
        config->signal,
        sizeof(float)*config->window_length);
    memcpy(
        ret->signal_end,
        config->signal+config->signal_length-config->window_length,
        sizeof(float)*config->window_length);
    ret->config = *config;
    return ret;
fail:
    windowed_lookup_f32_free(ret);
    return NULL;
}
    
void
windowed_lookup_f32_out_of_bounds_gen_const(float *dest, unsigned int dest_length, void *aux)
{
    float *val = aux;
    while (dest_length-- > 0) {
        dest[dest_length] = *val;
    }
}

struct windowed_lookup_f32_access_t windowed_lookup_f32_access (
    struct windowed_lookup_f32_t *wlu, int index)
{
    const float *signal_section;
    int W = wlu->config.window_length,
            signal_length = wlu->config.signal_length;
    if (index < -W) { index = -W; }
    if (index > signal_length) { index = signal_length; }
    if (index < 0) { signal_section = wlu->signal_start + (index + W); }
    else if (index >= (signal_length - W)) {
        signal_section = wlu->signal_end + (W - (signal_length - index));
    } else {
        signal_section = wlu->config.signal + index;
    }
    struct windowed_lookup_f32_access_t ret = {
        .signal_section = signal_section,
        .adjusted_index = index
    };
    return ret;
}
