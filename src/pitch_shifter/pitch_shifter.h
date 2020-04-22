#ifndef PITCH_SHIFTER_H
#define PITCH_SHIFTER_H 

#include <stdint.h>
#include "dsp_math.h"

struct pitch_shifter_config {
    const float *(*get_samples)(const u48q16 sample_index, void *aux);
    void *get_samples_aux;
    float ps_min;
    float ps_max;
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
};

struct pitch_shifter *
pitch_shifter_new(struct pitch_shifter_config *config);

void
pitch_shifter_process(struct pitch_shifter *self,
                      const u16q16 *ps_rate_sig,
                      const s16q16 *ts_rate_sig,
                      float *yi);

void pitch_shifter_set_position_at_block_start(
struct pitch_shifter *self,
u48q16 position);

void
pitch_shifter_clamp_ps_rate_sig(struct pitch_shifter *self,
u16q16 *ps_rate_sig, uint32_t length);

uint32_t
pitch_shifter_B(struct pitch_shifter *ps);

#endif /* PITCH_SHIFTER_H */
