#include <complex.h>
#include <math.h>
#include "dsp_math.h"
#include "kiss_fftr.h"

#define SWAP(a,b)\
    ({ typeof(a) tmp = a;\
       a = b;\
       b = tmp; })

void
dspm_mul_vf32_vf32_vf32(const float *src0,
                        const float *src1,
                        float *dst,
                        uint32_t length)
{
    while (length--) {
        *dst++ = *src0++ * *src1++;
    }
}

void
dspm_mul_vf32_f32(float *srcdst,
                  float src1,
                  uint32_t length)
{
    while (length--) {
        *srcdst++ *= src1;
    }
}

void
dspm_neg_vf32(float *srcdst,
              uint32_t length)
{
    /* Negation of IEEE floating point involves flipping highest-bit */
    while (length--) {
        *(uint32_t*)srcdst++ ^= 0x80000000;
    }
}

float
dspm_max_vf32(const float *src,
              uint32_t length)
{
    float ret = *src++;
    length--;
    while (length--) {
        ret = *src > ret ? *src : ret;
        src++;
    }
    return ret;
}

uint32_t
dspm_max_vu32(const uint32_t *src,
              uint32_t length)
{
    uint32_t ret = *src++;
    length--;
    while (length--) {
        ret = *src > ret ? *src : ret;
        src++;
    }
    return ret;
}

void
dspm_abs_vf32(float *srcdst, uint32_t length)
{
    while (length--) {
        *srcdst = fabsf(*srcdst);
        srcdst++;
    }
}

void
dspm_abs_vz32(float complex *srcdst, uint32_t length)
{
    while (length--) {
        *srcdst = cabsf(*srcdst);
        srcdst++;
    }
}

void
dspm_abs_vz32_vf32(const float complex *src, float *dst, uint32_t length)
{
    while (length--) {
        *dst++ = cabsf(*src++);
    }
}

void
dspm_sub_vf32_vf32_vf32(const float *src0,
                        const float *src1,
                        float *dst,
                        uint32_t length)
{
    while (length--) {
        *dst++ = *src0++ - *src1++;
    }
}


void
dspm_sub_vf32_vu32_vf32(const float *src0,
                        const uint32_t *src1,
                        float *dst,
                        uint32_t length)
{
    while (length--) {
        *dst++ = *src0++ - *src1++;
    }
}


void
dspm_add_vf32_f32_vf32(const float *src0, float src1, float *dst, uint32_t length)
{
    while (length--) {
        *dst++ = *src0++ + src1;
    }
}

void
dspm_add_vf32_vf32(float *srcdst, const float *src, uint32_t length)
{
    while (length--) {
        *srcdst++ += *src++;
    }
}

void
dspm_add_vu32_u32(uint32_t *srcdst,
                  uint32_t length,
                  uint32_t c)
{
    while (length--) {
        *srcdst++ += c;
    }
}

void
dspm_mul_vu32_u32(uint32_t *srcdst,
                  uint32_t length,
                  uint32_t c)
{
    while (length--) {
        *srcdst++ *= c;
    }
}

void
dspm_clip_below_vf32_f32(float *srcdst,
                         float lb,
                         uint32_t length)
{
    while (length--) {
        *srcdst = *srcdst < lb ? lb : *srcdst;
        srcdst++;
    }
}

/* Note that this is not compensated summation */
float
dspm_sum_vf32(const float *src, uint32_t length)
{
    float ret = 0;
    while (length--) {
        ret += *src++;
    }
    return ret;
}

/* Assumes d is non-zero */
void
dspm_div_vf32_f32(float *srcdst,
                  float d,
                  uint32_t length)
{
    while (length--) {
        *srcdst++ /= d;
    }
}

struct dspm_rfft_vf32_vz32_cfg {
    kiss_fftr_cfg cfg;
    uint32_t inverse;
};

void
dspm_rfft_vf32_vz32_cfg_free(struct dspm_rfft_vf32_vz32_cfg *cfg)
{
    if (cfg) {
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
                    float complex *freq)
{
    if (cfg->inverse) {
        kiss_fftri(cfg->cfg,freq,time);
    } else {
        kiss_fftr(cfg->cfg,time,freq);
    }
}

float
dspm_mean_vf32_f32(const float *src, uint32_t length)
{
    float ret = dspm_sum_vf32(src, length);
    ret /= length;
    return ret;
}

void
dspm_rev_vu32(uint32_t *srcdst, uint32_t length)
{
    uint32_t beg = 0, end = length - 1;
    while (beg < (length>>1)) {
        SWAP(srcdst[beg],srcdst[end]);
        beg++;
        end--;
    }
}

void
dspm_floor_vf32_vu32(const float *src, uint32_t *dst, uint32_t length)
{
    while (length--) {
        *dst++ = (uint32_t)floor(*src++);
    }
}

void
dspm_lookup_vf32_vu32_vf32(const float *src0, const uint32_t *src1,
float *dst, uint32_t length)
{
    while (length--) {
        *dst++ = src0[*src1++];
    }
}

void
dspm_floor_vu24q8_vu32(const u24q8 *src,
                       uint32_t *dst,
                       uint32_t N)
{
    while (N--) {
        *dst++ = *src++ >> 8;
    }
}

/* note: result will be wrong if src1 contains values >= 2^24 */
void
dspm_sub_vu24q8_vu32_vf32(const u24q8 *src0,
                          const uint32_t *src1,
                          float *dst,
                          uint32_t N)
{
    const float scale = 1./256.;
    while (N--) {
        *dst++ = (*src0++ - (*src1++ << 8))*scale;
    }
}

static inline void
cumsum_vu32_u64_vu64(const uint32_t *src,
                     uint64_t initial_sum,
                     uint64_t *dst,
                     uint32_t N)
{
    while (N--) {
        initial_sum += *src++;
        *dst++ = initial_sum;
    }
}

void
dspm_cumsum_vu16q16_u48q16_vu48q16(const u16q16 *src,
                                   u48q16 initial_sum,
                                   u48q16 *dst,
                                   uint32_t N)
{
    cumsum_vu32_u64_vu64(src,initial_sum,dst,N);
}

struct dspm_2dline_s48q16
dspm_2dline_s48q16_points(s48q16 x0,
        s48q16 y0,
        s48q16 x1,
        s48q16 y1)
{
    struct dspm_2dline_s48q16 ret = {
        .x0 = x0;
        .m = ((y2-y1)<<16)/(x1-x0);
        .b = y1;
    };
}

void
dspm_2dline_s48q16_lookup_vs48q16(const struct dspm_2dline_s48q16 *restrict line,
                                  s48q16 *restrict x,
                                  uint32_t N)
{
    while (N--) {
        /* y = m*(x - x0) + y */
        *x = ((line->m * (x - line->x0)) >> 16) + line->y;
        x++;
    }
}

