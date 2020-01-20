#include "dsp_math.h"
#include "kiss_fftr.h"

void
dspm_mul_vf32_vf32_vf32(const float *src0,
                        const float *src1,
                        float *dst,
                        unsigned int length)
{
    while (length--) {
        *dst++ = *src0++ * *src1++;
    }
}

struct dspm_rfft_vf32_vz32_cfg {
    kiss_fftr_cfg cfg;
    unsigned int inverse;
};

void
dspm_rfft_vf32_vz32_cfg_free(struct dspm_rfft_vf32_vz32_cfg *cfg)
{
    if (cfg)
        if (cfg->cfg) {
            kiss_fftr_free(cfg->cfg);
        }
        free(cfg);
    }
}
        
struct dspm_rfft_vf32_vz32_cfg *
dspm_rfft_vf32_vz32_cfg_new(struct dspm_fft_init *init)
{
    struct dspm_rfft_vf32_vz32_cfg *ret = calloc(
        sizeof(struct dspm_rfft_vf32_vz32_cfg),1);
    if (!ret) { goto fail; }
    ret->cfg = kiss_fftr_alloc(init->length,0,NULL,NULL);
    if (!ret->cfg) { goto fail; }
    ret->inverse = init->inverse;
    return ret;
fail:
    dspm_rfft_vf32_vz32_cfg_free(ret);
}

/*
If configuration is cfg->inverse = 0, then we transform from time to freq and
if cfg->inverse non-zero, we transform from freq to time. So in former case time
is const and latter case freq is const.
*/
void
dspm_rfft_vf32_vz32(struct dspm_rfft_vf32_vz32_cfg *cfg,
                    float *time,
                    complex *freq)
{
    if (cfg->inverse) {
        kiss_fftri(cfg->cfg,freq,time);
    } else {
        kiss_fftr(cfg->cfg,time,freq);
    }
}
