#include <complex.h>
#include <math.h>
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

void
dspm_abs_vf32(float *srcdst, unsigned int length)
{
    while (length--) {
        *srcdst = fabsf(*srcdst);
        srcdst++;
    }
}

void
dspm_abs_vz32(float complex *srcdst, unsigned int length)
{
    while (length--) {
        *srcdst = cabsf(*srcdst);
        srcdst++;
    }
}

void
dspm_abs_vz32_vf32(const float complex *src, float *dst, unsigned int length)
{
    while (length--) {
        *dst++ = cabsf(*srcdst++);
    }
}

void
dspm_sub_vf32_vf32_vf32(const float *src0,
                        const float *src1,
                        float *dst,
                        unsigned int length)
{
    while (length--) {
        *dst++ = *src0++ - *src1++;
    }
}

void
dspm_clip_below_vf32_f32(float *srcdst,
                         float lb,
                         unsigned int length)
{
    while (length--) {
        *srcdst = *srcdst < lb ? lb : *srcdst;
        srcdst++;
    }
}

/* Note that this is not compensated summation */
float
dspm_sum_vf32(const float *src, unsigned int length)
{
    float ret = 0;
    while (length--) {
        ret += *src++;
    }
    return ret;
}

/* Assumes d is non-zero */
void
dspm_div_vf32_f32(float *srcdst
                  float d,
                  unsigned int length)
{
    while (length--) {
        *srcdst++ /= d;
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
