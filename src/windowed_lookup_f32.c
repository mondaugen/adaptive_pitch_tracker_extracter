/*
Data structure for lookups of f32 that will have a fixed window size and that
may request out of bounds data.
*/

#include <string.h>

struct windowed_lookup_init_t {
    /* The size of the lookups */
    unsigned int window_length;
    /* The signal to lookup from. This is not freed when windowed_lookup_t is freed. */
    const float *signal;
    /* The length of this signal */
    unsigned int signal_length;
    /* A function that allows you to fill the out-of-bounds values. This
    function is only called on initialization. If this is NULL, the
    out-of-bounds values are filled with 0s. */
    void (*out_of_bounds_gen)(
        /* Where the values will go */
        float *dest,
        /* How many values to fill */
        unsigned int dest_length,
        /* Auxilary data */
        void *aux);
    /* Auxilary data for out_of_bounds_gen. This is not freed when
    windowed_lookup_t is freed. */    
    void *aux;
};

struct windowed_lookup_t {
    struct windowed_lookup_init_t config;
    /* An array of length 2 * signal_length holding signal_length fill values
    followed by the first signal_length values */
    float *signal_start;
    /* An array of length 2 * signal_length holding the last signal_length
    values followed by signal_length fill values */
    float *signal_end;
};

struct windowed_lookup_access_t {
    /* Section of signal of length window_length */
    const float *signal_section;
    /* Index after access function adjusted it, in the case that the access was
    out of bounds. So the minimum value this could take is -window_length and
    the maxmimum signal_length. */
    const int adjusted_index;
};

static int init_arg_chk(struct windowed_lookup_init_t *wli)
{
    if (!wli->signal) { return -1 }
    if (wli->window_length < 0) { return -2; }
    if (wli->signal_length < 0) { return -3; }
    return 0;
}

void
windowed_lookup_free(struct windowed_lookup_t *wl)
{
    if (!wl) { return; }
    if (wl->signal_start) { free(wl->signal_start); }
    if (wl->signal_end) { free(wl->signal_end); }
    free(wl);
}

struct windowed_lookup_t *
windowed_lookup_new(windowed_lookup_init_t *config)
{
    if (init_arg_chk(config)) { return NULL; }
    windowed_lookup_t *ret = calloc(1,sizeof(struct windowed_lookup_t));
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
    windowed_lookup_free(ret);
    return NULL;
}
    
void
windowed_lookup_out_of_bounds_gen_const(float *dest, unsigned int dest_length, void *aux)
{
    float *val = aux;
    while (dest_length-- > 0) {
        dest[dest_length] = *val;
    }
}

void
windowed_lookup_access (struct windowed_lookup_t *wlu, int index, struct windowed_lookup_access_t *wla)
{
    unsigned int W = wlu->config.window_length,
                 signal_length = wlu->config.signal_length;
    if (index < -W) { index = -W; }
    if (index > signal_length) { index = signal_length; }
    if (index < 0) { wla->signal_section = wlu->signal_start + (index + W); }
    else if (index >= (signal_length - W)) {
        wla->signal_section = wlu->signal_start + (W - (signal_length - index));
    } else {
        wla->signal_section = wlu->config.signal + index;
    }
    wla->adjusted_index = index;
}
