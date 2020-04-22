/* Return samples of a sine wave */
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include "pitch_shifter.h"
#include "dsp_math.h"

struct ps_sine_samples {
    /* normalized frequency */
    float f;
    /* block size */
    uint32_t B;
    /* phase */
    float p;
    float *samples;
};

static void
ps_sine_samples_free(struct ps_sine_samples *ret)
{
    if (ret) {
        if (ret->samples) { free(ret->samples); }
        free(ret);
    }
}

static struct ps_sine_samples *
ps_sine_samples_new(float f, uint32_t B)
{
    struct ps_sine_samples *ret = calloc(1,sizeof(struct ps_sine_samples));
    if (!ret) { return NULL; }
    ret->samples = calloc(B,sizeof(float));
    if (!ret->samples) { goto fail; }
    ret->B = B;
    ret->f = f;
    ret->p = 0;
    return ret;
fail:
    ps_sine_samples_free(ret);
    return NULL;
}

static const float *
get_samples(const u48q16 sample_index, void *aux)
{
    /* we have no signal to index so we just ignore sample_index */
    struct ps_sine_samples *self = aux;
    uint32_t n;
    for (n = 0; n < self->B; n++) {
        self->samples[n] = sin(2.*M_PI*self->p);
        self->p += self->f;
    }
    while (self->p >= 1.) { self->p -= 1.; }
    while (self->p < 0.) { self->p += 1.; }
    return self->samples;
}

/*
f is normalized frequency, initialize B in config before calling this, returns
non-zero on error
*/
int
ps_sine_samples_config(float f, struct pitch_shifter_config *config)
{
    config->get_samples_aux = ps_sine_samples_new(f,config->B);
    if (!config->get_samples_aux) {
        return -1;
    }
    config->get_samples = get_samples;
    return 0;
}
