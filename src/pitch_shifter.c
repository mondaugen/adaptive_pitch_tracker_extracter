#include <stdint.h>
#include <stdlib.h>
#include "rngbuf_f32.h"
#include "dsp_math.h"

struct pitch_shifter {
    struct rngbuf_f32 *sig_rb;
    int sig_rb_idcs_valid;
    int64_t sig_rb_min_idx;
    int64_t sig_rb_max_idx;
    const float *(*get_samples)(const u48q16 sample_index, void *aux);
    void *get_samples_aux;
    float ps_min;
    float ps_max;
    uint32_t B;
    void (*interpolator)(const u24q8 *xi,
                         const float *y,
                         float *yi,
                         uint32_t N,
                         void *aux);
    void (*interpolator_range)(const u48q16 pos_first,
                               const u48q16 pos_last,
                               int64_t* req_first,
                               int64_t* req_last,
                               void* aux);
    uint32_t (*interpolator_n_points)(uint32_t N, void* aux);
    void *interpolator_aux;
    s48q16 time_at_block_start;
    s48q16 pos_at_block_start;
};

static int check_pitch_shifter_config(struct pitch_shifter_config *config)
{
    /*
    We restrict ps_max to 12 octaves.
    B must at least be limited by the (2^32 - 1)/config->ps_max because that
    is the maxmimum size of a ring buffer and we look up at most
    B*config->ps_max per processing block, so the most that could be discarded
    is B*config->ps_max. Therefore limiting B to 2^16 should be safe.
    */
    if (!config) { return -1; }
    if (!config->get_samples) { return -2; }
    if (config->ps_min <= 0) { return -3; }
    if ((config->ps_max < config->ps_min)
        || (config->ps_max <= 0)
        || (config->ps_max > (1 << 12)) { return -4; }
    if ((config->B > (1 << 16)) || (config->B <= 0)) { return -5; }
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
    uint32_t sig_rb_size = 0;
    struct pitch_shifter *self = calloc(1,sizeof(struct pitch_shifter));
    if (!self) { goto fail; }
    // The most values that could be requested are
    // ceil((B*ps_max+get_interpolator_n_points(B))/B)*B
    tmp=config->B*config->ps_max+
        config->get_interpolator_n_points(B,config->interpolator_aux);
    while sig_rb_size < tmp:
        sig_rb_size += config->B
    self->sig_rb = rngbuf_f32_new(sig_rb_size);
    if (!self->sig_rb) { goto fail; }
    // These values are invalid until after the first call to process
    // This is the index of the value at index 0 in the ring buffer

    self->sig_rb_idcs_valid = 0;
    self->sig_rb_min_idx = -1;
    // This is the index of the value at index
    // self.sig_rb.contents_size() - 1 in the ringbuffer
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

void pitch_shifter_set_position_at_block_start(u48q16 position)
{
    self->pos_at_block_start = position;
    self->time_at_block_start = position;
    self->sig_rb_idcs_valid = 0;
}

/*
ps_pos_sig and ts_pos_sig have length self->B+1
enough values are fetched using get_samples and the times in ts_pos_sig,
the interpolated values are returned in y (which has length self->B) to obtain the
pitch shift.
*/
static void
process_pos_sig(struct pitch_shifter *self,
                const u48q16 *ps_pos_sig,
                const  u48q16 *ts_pos_sig,
                float *y)
{
    int64_t first_required_idx, last_required_idx, fetch_start_idx;
    uint32_t n_discard;
    self->interpolator_range(ps_pos_sig[0],ps_pos_sig[-2],
    &first_required_idx,&last_required_idx);

    if (self->sig_rb_idcs_valid) {
        /* self->sig_rb_max_idx's minimum value is -1 */
        fetch_start_idx = self->sig_rb_max_idx + 1;
        /* discard all values between self.sig_rb_min_idx and the first_required_idx */
        n_discard = first_required_idx - self->sig_rb_min_idx;
        rngbuf_f32_advance_head(self->sig_rb,n_discard);
    } else {
        rngbuf_reset((struct rngbuf *)self->sig_rb);
        fetch_start_idx = first_required_idx;
        /* start here so the accumlated value is correct (see the while loop below) */
        self->sig_rb_max_idx = fetch_start_idx - 1;
        self->sig_rb_idcs_valid = 1;
    }
    self.sig_rb_min_idx = first_required_idx

    n_fetch_vals = last_required_idx - fetch_start_idx + 1

    /* linear interpolation giving the look up times from the signal indices */
    fetch_time_interpolator=interpolate.interp1d(
        [ps_pos_sig[0],ps_pos_sig[-1]],
        [ts_pos_sig[0],ts_pos_sig[-1]],
        fill_value='extrapolate')

}
