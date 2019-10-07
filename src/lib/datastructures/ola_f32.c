/*
Overlap and add buffer for 32bit floats.
*/

#include <stdlib.h>
#include "ola_f32.h"

struct ola_f32_t {
    struct ola_f32_init_t config;
    float *buffer;
    unsigned int len_mask;
    unsigned int offset;
};

/* Sum in sum_in_length values into a pvs_ola_t */
void
ola_f32_sum_in(struct ola_f32_t *ola, const float *input)
{
    ola_f32_add(ola->buffer + ola->offset, input, ola->config.sum_in_length - ola->offset);
    ola_f32_add(ola->buffer, input + ola->config.sum_in_length - ola->offset, ola->offset);
}

/* Return pointer to contigious region holding shift_out_length values */
const float *
ola_f32_shift_out(struct ola_f32_t *ola)
{
    const float *ret = ola->buffer + ola->offset;
    ola->offset = (ola->offset + ola->config.shift_out_length) & ola->len_mask;
    return ret;
}

static int
chk_pow_2(unsigned int x)
{
    unsigned int n = 1;
    while (n < x) { n <<= 1; }
    if (n == x) { return 1; }
    return 0;
}
    
static int
config_chk(struct ola_f32_init_t *config)
{
    if (!config) { return -1; }
    /* sum_in_length must be power of 2 */
    if (!chk_pow_2(config->sum_in_length)) { return -2; }
    /* the shift_out_length must be divisor of sum_in_length */
    if ((config->sum_in_length % config->shift_out_length) != 0) { return -3; }
    return 0;
}

void
ola_f32_free(struct ola_f32_t *ola)
{
    if (!ola) { return; }
    if (ola->buffer) { free(ola->buffer); }
    free(ola);
}

struct ola_f32_t *
ola_f32_new(struct ola_f32_init_t *config)
{ 
    int config_chk_ret = config_chk(config);
    if (config_chk_ret) { return NULL; }
    struct ola_f32_t *ret = calloc(1,sizeof(struct ola_f32_t));
    if (!ret) { goto fail; }
    ret->config = *config;
    ret->buffer = calloc(1,sizeof(float)*ret->config.sum_in_length);
    if (!ret->buffer) { goto fail; }
    ret->len_mask = ret->config.sum_in_length - 1;
    /* offset is already 0 */
    return ret;
fail:
    ola_f32_free(ret);
    return NULL;
}
