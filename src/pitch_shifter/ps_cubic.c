/*
Get a configuration that does cubic interpolation when pitch shifting.
*/

#include "ps_cubic.h"
#include "dsp_math.h"

static void 
interpolator(const u16q16 *xi,
                  const float *y,
                  float *yi,
                  uint32_t N,
                  void *aux)
{
    dspm_interp1d4p_vu16q16_vf32_vf32(xi,y,yi,N);
}

static void
interpolator_range(u48q16 pos_first,
                               u48q16 pos_last,
                               int64_t* req_first,
                               int64_t* req_last,
                               void* aux)
{
    /* def default_get_interpolator_range(x_min,x_max):
       return (int(np.floor(x_min)-1),int(np.floor(x_max)+2)) */
    *req_first = (pos_first >> 16) - 1;
    *req_last = (pos_last >> 16) + 2;
}

static uint32_t
interpolator_n_points(uint32_t N, void* aux)
{
    return N+3;
}

/* prefill interpolation fields of configuration structure */
void
ps_cubic_config(struct pitch_shifter_config *config)
{
    config->interpolator = interpolator;
    config->interpolator_range = interpolator_range;
    config->interpolator_n_points = interpolator_n_points;
}
