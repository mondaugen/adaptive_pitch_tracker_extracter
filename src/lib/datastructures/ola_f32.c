/*
Overlap and add buffer for 32bit floats.
*/

#include <stdlib.h>
#include <string.h>
#include "ola_f32.h"

struct ola_f32_t {
    struct ola_f32_init_t config;
    float *buffer;
    unsigned int buffer_len;
    unsigned int len_mask;
    unsigned int offset;
};

/*
Sum in sum_in_length values into a pvs_ola_t and return pointer to region of
length shift_out_length containing the most recently complete overlap-and-add
frame (the frame into which no more additions will occur). If this frame is not
used before the next call to ola_f32_sum_in_and_shift_out, it will be set to 0
the next call to ola_f32_sum_in_and_shift_out.
*/
const float *
ola_f32_sum_in_and_shift_out(struct ola_f32_t *ola, const float *input)
{
    /* Zero last frame */
    unsigned int zero_start_0 = (ola->offset - ola->config.shift_out_length) & ola->len_mask,
                 /* If zero_start_0 comes after the current offset, it means
                 zero_start_0 was wrapped around, so its first region is valid
                 from zero_start_0 to the end of the buffer. Otherwise it's just
                 the shift_out_length. */
                 zero_len_0 = (zero_start_0 > ola->offset) ? ola->buffer_len - zero_start_0 :
                    ola->config.shift_out_length,
                 zero_len_1 = ola->config.shift_out_length - zero_len_0;
    memset(ola->buffer + zero_start_0, 0, sizeof(float) * zero_len_0);
    /* NOTE: memset should not fail, but just do nothing, with length of 0 */
    memset(ola->buffer, 0, sizeof(float) * zero_len_1);
    /* Now sum in the input */
    unsigned int sum_len_0 = 
        ((ola->offset + ola->config.sum_in_length) & ola->len_mask) < ola->offset ?
        ola->buffer_len - ola->offset : ola->config.sum_in_length,
                 sum_len_1 = ola->config.sum_in_length - sum_len_0;
    ola_f32_add(ola->buffer + ola->offset, input, sum_len_0);
    /* ola_f32_add should not fail, but just do nothing, with length of 0 */
    ola_f32_add(ola->buffer, input + sum_len_0, sum_len_1);
    /* Return start of the most recently finished region */
    const float *ret = ola->buffer + ola->offset;
    ola->offset = (ola->offset + ola->config.shift_out_length) & ola->len_mask;
    return ret;
}

static int
next_pow_2(unsigned int x)
{
    unsigned int n = 1;
    while (n < x) { n <<= 1; }
    return n;
}
    
static int
config_chk(struct ola_f32_init_t *config)
{
    if (!config) { return -1; }
    if (config->sum_in_length < 1) { return -2; }
    if (config->shift_out_length > config->sum_in_length) { return -3; }
    /* We always want to return contiguous memory, so we ensure the
    shift_out_length divides the buffer size */
    if (next_pow_2(config->shift_out_length + config->sum_in_length) 
        % config->shift_out_length) { return -4; }
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
    ret->buffer_len = next_pow_2(
        ret->config.sum_in_length + ret->config.shift_out_length);
    ret->buffer = calloc(ret->buffer_len,sizeof(float));
    if (!ret->buffer) { goto fail; }
    ret->len_mask = ret->buffer_len - 1;
    /* offset is already 0 */
    return ret;
fail:
    ola_f32_free(ret);
    return NULL;
}
