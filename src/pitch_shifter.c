#include <stdlib.h>
#include "rngbuf_f32.h"
#include "dsp_math.h"

struct pitch_shifter {
    struct rngbuf_f32 *sig_rb;
    int sig_rb_idcs_valid;
    int sig_rb_min_idx;
    int sig_rb_max_idx;
    const float *(*get_samples)(const u24q8 sample_index, void *aux);
    void *get_samples_aux;
    float ps_min;
    float ps_max;
    unsigned int B;
    void (*interpolator)(const u24q8 *xi,
                         const float *y,
                         float *yi,
                         unsigned int N,
                         void *aux);
    void (*interpolator_range)(const u24q8 pos_first,
                               const u24q8 pos_last,
                               int* req_first,
                               int* req_last,
                               void* aux);
    unsigned int (*interpolator_n_points)(unsigned int N, void* aux);
    void *interpolator_aux;
    u24q8 time_at_block_start;
    u24q8 pos_at_block_start;
};

static int check_pitch_shifter_config(struct pitch_shifter_config *config)
{
    if (!config) { return -1; }
    if (!config->get_samples) { return -2; }
    if (config->ps_min <= 0) { return -3; }
    if ((config->ps_max < config->ps_min) || (config->ps_max <= 0)) { return -4; }
    if (config->B <= 0) { return -5; }
    if (!config->interpolator) { return -6; }
    if (!config->interpolator_range) { return -7; }
    if (!config->interpolator_n_points) { return -8; }
    return 0;
}

struct pitch_shifter *
pitch_shifter_new(struct pitch_shifter_config *config)
{
    int err = check_pitch_shifter_config(config);
    if (err) { return NULL; }
    float tmp;
    unsigned int sig_rb_size = 0;
    struct pitch_shifter *self = calloc(1,sizeof(struct pitch_shifter));
    if (!self) { goto fail; }
    //    # The most values that could be requested are
    //    # ceil((B*ps_max+get_interpolator_n_points(B))/B)*B
    tmp=config->B*config->ps_max+
        config->get_interpolator_n_points(B,config->interpolator_aux);
    while sig_rb_size < tmp:
        sig_rb_size += config->B
    self->sig_rb = rngbuf_f32_new(sig_rb_size);
    if (!self->sig_rb) { goto fail; }
    //    # These values are invalid until after the first call to process
    //    # This is the index of the value at index 0 in the ring buffer

    self->sig_rb_idcs_valid = 0;
    self->sig_rb_min_idx = -1;
    //    # This is the index of the value at index
    //    # self.sig_rb.contents_size() - 1 in the ringbuffer
    self->sig_rb_max_idx = -1;
    self->get_samples = get_samples;
    self->get_samples_aux = get_samples_aux;
    self->ps_min = ps_min;
    self->ps_max = ps_max;
    self->B = B;
    self->interpolator=interpolator;
    self->interpolator_range=interpolator_range;
    self->interpolator_aux=interpolator_aux;
    self->time_at_block_start=0;
    self->pos_at_block_start=0;
fail:
    if (self) {
        if (self->sig_rb) { rngbuf_f32_free(self->sig_rb); }
        free(self);
    }
    return NULL;
}

void pitch_shifter_set_position_at_block_start(u24q8 position)
{
    self->pos_at_block_start = position;
    self->time_at_block_start = position;
    self->sig_rb_idcs_valid = 0;
}
