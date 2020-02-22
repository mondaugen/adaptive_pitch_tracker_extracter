#include "dsp_math.h"

void
spec_diff_dft_free(void *aux)
{
    dspm_rfft_vf32_vz32_cfg_free(aux);
}

void *
spec_diff_dft_new(unsigned int W)
{
    struct dspm_fft_init init = {
        .length = W,
        .inverse = 0
    };
    return dspm_rfft_vf32_vz32_cfg_new(&init);
}

void spec_diff_forward_dft(void *aux, float *time, float complex *freq)
{
    struct dspm_rfft_vf32_vz32_cfg *cfg = aux;
    dspm_rfft_vf32_vz32(cfg,time,freq);
}

unsigned int spec_diff_complex_size(unsigned int W)
{
    /* TODO: shouldn't dsp_math have a routine to get this value? */
    return W/2 + 1;
}
