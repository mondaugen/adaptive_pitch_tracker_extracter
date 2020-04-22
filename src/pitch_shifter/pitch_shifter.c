#include <stdint.h>
#include <stdlib.h>
#include "pitch_shifter.h"
#include "ringbuffer.h"
#include "rngbuf_f32.h"

struct pitch_shifter {
    struct rngbuf_f32 *sig_rb;
    int sig_rb_idcs_valid;
    int64_t sig_rb_min_idx;
    int64_t sig_rb_max_idx;
    const float *(*get_samples)(const u48q16 sample_index, void *aux);
    void *get_samples_aux;
    u16q16 ps_min;
    u16q16 ps_max;
    uint32_t B;
    void (*interpolator)(const u16q16 *xi,
                         const float *y,
                         float *yi,
                         uint32_t N,
                         void *aux);
    void (*interpolator_range)(u48q16 pos_first,
                               u48q16 pos_last,
                               int64_t* req_first,
                               int64_t* req_last,
                               void* aux);
    uint32_t (*interpolator_n_points)(uint32_t N, void* aux);
    void *interpolator_aux;
    s48q16 time_at_block_start;
    u48q16 pos_at_block_start;
};

uint32_t
pitch_shifter_B(struct pitch_shifter *ps) { return ps->B; }

static int check_pitch_shifter_config(struct pitch_shifter_config *config)
{
    /*
    We restrict ps_max to 3 octaves.
    B*ps_max must at least be limited by 2**16 because we use u16q16 as a
    relative position index and we look up at most B*config->ps_max per
    processing block, so the most that could be discarded is B*config->ps_max.
    Therefore we limit B to 2**16/ps_max. Don't forget the interpolator might
    want some extra values, so that's why we limit to 3 octaves, not 4.
    */
    if (!config) { return -1; }
    if (!config->get_samples) { return -2; }
    if (config->ps_min <= 0) { return -3; }
    if ((config->ps_max < config->ps_min)
        || (config->ps_max <= 0)
        || (config->ps_max > (1 << 3))) { return -4; }
    if (!config->interpolator) { return -6; }
    if (!config->interpolator_range) { return -7; }
    if (!config->interpolator_n_points) { return -8; }
    if ((config->interpolator_n_points(config->B*config->ps_max,
         config->interpolator_aux) > (1 << 16)) || (config->B <= 0)) {
         return -5;
    }
    return 0;
}

static inline u16q16
float_to_u16q16(float f)
{
    u16q16 ret;
    dspm_cvt_vf32_vu16q16(&f, &ret, 1);
    return ret;
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
        config->interpolator_n_points(config->B,config->interpolator_aux);
    while (sig_rb_size < tmp) {
        sig_rb_size += config->B;
    }
    self->sig_rb = rngbuf_f32_new(sig_rb_size);
    if (!self->sig_rb) { goto fail; }
    // These values are invalid until after the first call to process
    // This is the index of the value at index 0 in the ring buffer

    self->sig_rb_idcs_valid = 0;
    self->sig_rb_min_idx = -1;
    // This is the index of the value at index
    // self.sig_rb.contents_size() - 1 in the ringbuffer
    self->sig_rb_max_idx = -1;
    self->get_samples = config->get_samples;
    self->get_samples_aux = config->get_samples_aux;
    self->ps_min = float_to_u16q16(config->ps_min);
    self->ps_max = float_to_u16q16(config->ps_max);
    self->B = config->B;
    self->interpolator=config->interpolator;
    self->interpolator_range=config->interpolator_range;
    self->interpolator_aux=config->interpolator_aux;
    self->time_at_block_start=0;
    self->pos_at_block_start=0;
    return self;
fail:
    if (self) {
        if (self->sig_rb) { rngbuf_f32_free(self->sig_rb); }
        free(self);
    }
    return NULL;
}

void pitch_shifter_set_position_at_block_start(struct pitch_shifter *self, u48q16 position)
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
                u48q16 *ps_pos_sig,
                s48q16 *ts_pos_sig,
                float *yi)
{
    int64_t first_required_idx, last_required_idx, fetch_start_idx;
    uint32_t n_discard, n_fetch_vals, n;
    s48q16 fetch_time;
    u16q16 local_ps_pos_sig[self->B];
    self->interpolator_range(ps_pos_sig[0],ps_pos_sig[-2],
    &first_required_idx,&last_required_idx,self->interpolator_aux);

    if (self->sig_rb_idcs_valid) {
        /* self->sig_rb_max_idx's minimum value is -1 */
        fetch_start_idx = self->sig_rb_max_idx + 1;
        /* discard all values between self.sig_rb_min_idx and the first_required_idx */
        n_discard = first_required_idx - self->sig_rb_min_idx;
        rngbuf_f32_advance_head(self->sig_rb,n_discard);
    } else {
        rngbuf_reset((struct rngbuf *)self->sig_rb);
        fetch_start_idx = first_required_idx;
        /* start here so the accumulated value is correct (see the while loop below) */
        self->sig_rb_max_idx = fetch_start_idx - 1;
        self->sig_rb_idcs_valid = 1;
    }
    self->sig_rb_min_idx = first_required_idx;

    n_fetch_vals = last_required_idx - first_required_idx + 1;

    /* linear interpolation giving the look up times from the signal indices */
    struct dspm_2dline_s48q16
    fetch_time_interpolator=dspm_2dline_s48q16_points(
    ps_pos_sig[0],ts_pos_sig[0],ps_pos_sig[self->B],ts_pos_sig[self->B]);
    
    while (fetch_start_idx < last_required_idx) { 
        fetch_time = fetch_start_idx << 16;
        dspm_2dline_s48q16_lookup_vs48q16(&fetch_time_interpolator,
                                          &fetch_time,
                                          1);
        const float *samples = self->get_samples(fetch_time,self->get_samples_aux);
        rngbuf_f32_push_copy(self->sig_rb, samples, self->B);
        fetch_start_idx += self->B;
        self->sig_rb_max_idx += self->B;
    }
    /* offset by first_required_idx to get a local index signal */
    dspm_sub_vu48q16_s64_vu16q16(ps_pos_sig,
                                 first_required_idx,
                                 local_ps_pos_sig,
                                 self->B);
    /* get the required values */
    float y[n_fetch_vals];
    rngbuf_f32_memcpy(
        self->sig_rb,
        0,
        last_required_idx-first_required_idx+1,
        y);
    /* interpolator assumes domain is 0 to N - 1 */
    self->interpolator(local_ps_pos_sig,
                       y,
                       yi,
                       self->B,
                       self->interpolator_aux);
    self->time_at_block_start = ts_pos_sig[self->B];
    self->pos_at_block_start = ps_pos_sig[self->B];
}

void
pitch_shifter_clamp_ps_rate_sig(struct pitch_shifter *self,
u16q16 *ps_rate_sig, uint32_t length)
{
    /* limit pitch shift amount */
    dspm_min_vu16q16_u16q16(ps_rate_sig,self->ps_max,length);
    dspm_max_vu16q16_u16q16(ps_rate_sig,self->ps_min,length);
}

/* 
returns self->B samples of the signal (obtained through get_samples)
pitch-shifted according to ps_rate_sig and time-stretched according to
ts_rate_sig.
puts result in yi
yi, ps_rate_sig and ts_rate_sig must have length self->B
ps_rate_sig must be between self->ps_min and self->ps_max (.e.g, use
pitch_shifter_clamp_ps_rate_sig)
*/
void
pitch_shifter_process(struct pitch_shifter *self,
                      const u16q16 *ps_rate_sig,
                      const s16q16 *ts_rate_sig,
                      float *yi)
{
    u48q16 ps_pos_sig[self->B+1];
    s48q16 ts_pos_sig[self->B+1];
    ps_pos_sig[0] = self->pos_at_block_start;
    dspm_cumsum_vu16q16_u48q16_vu48q16(ps_rate_sig,
                                       self->pos_at_block_start,
                                       &ps_pos_sig[1],
                                       self->B);
    ts_pos_sig[0] = self->time_at_block_start;
    dspm_cumsum_vs16q16_s48q16_vs48q16(ts_rate_sig,
                                       self->pos_at_block_start,
                                       &ts_pos_sig[1],
                                       self->B);
    process_pos_sig(self,ps_pos_sig,ts_pos_sig,yi);
}
