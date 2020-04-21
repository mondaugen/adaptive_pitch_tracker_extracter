#ifndef DSP_MATH_H
#define DSP_MATH_H 

#include <stdint.h>
#include <complex.h>

/*
DSP math routines
The naming convention here is similar to Intel's IPP.
functions begin with dspm_
then a short description of the operation (e.g., dspm_add)
then the types, prepended optionally with v if a vector is accepted as an
argument (a pointer), followed by the type name (e.g., f32 for 32-bit float)
so the function dspm_add_vf32_vf32(float *srcdst, const float* src, uint32_t
length) computes srcdst[n] += src[n] for n in [0,length) and the function
dspm_add_vf32_vf32_vf32(const float *src0, const float *src1, float *dst,
uint32_t length) computs dst[n] = src0[n] + src1[n] for n in [0,length)

For operations requiring vectors, calling these on vectors of length zero or
NULL pointer arguments, the results are undefined.

*/

/* fixed-point types */
/* undefined with 16 integer and 16 fractional bits */
typedef uint32_t u16q16;
/* unsigned with 24 integer and 8 fractional bits */
typedef uint32_t u24q8;
/* signed with 24 integer and 8 fractional bits */
typedef int32_t s24q8;
/* unsigned with 48 integer and 16 fractional bits */
typedef uint64_t u48q16;
/* signed with 48 integer and 16 fractional bits */
typedef int64_t s48q16;

/* Fast floor(log2(x)) */
static inline int32_t
dspm_fast_floor_log2_f32(float x)
{
    uint32_t *ux = (uint32_t*)&x;
    /* Extracts exponent from IEEE 754 floating point number */
    int32_t exp = ((*ux>>23)&0xff) - 127;
    return exp;
}

/* Fast approximate fractional part of log2(x) */
static inline float
dspm_fast_log2_aprox_frac_f32(float x)
{
    uint32_t *ux = (uint32_t*)&x;
    /* (*ux&0x7fffff)**(2**-23) */
    float frac = (*ux&0x7fffff)*1.1920928955078125e-07;
    return frac;
}
    
void
dspm_mul_vf32_vf32_vf32(const float *src0,
                        const float *src1,
                        float *dst,
                        uint32_t length);

void
dspm_mul_vf32_f32(float *srcdst,
                  float src1,
                  uint32_t length);

void
dspm_neg_vf32(float *srcdst,
              uint32_t length);

float
dspm_max_vf32(const float *src,
              uint32_t length);

uint32_t
dspm_max_vu32(const uint32_t *src,
              uint32_t length);

void
dspm_add_vu32_u32(uint32_t *srcdst,
                  uint32_t length,
                  uint32_t c);

void
dspm_add_vf32_f32_vf32(const float *src0, float src1, float *dst, uint32_t length);

void
dspm_add_vf32_vf32(float *srcdst, const float *src, uint32_t length);

void
dspm_mul_vu32_u32(uint32_t *srcdst,
                  uint32_t length,
                  uint32_t c);

void
dspm_abs_vf32(float *srcdst, uint32_t length);

void
dspm_abs_vz32(float complex *srcdst, uint32_t length);

void
dspm_abs_vz32_vf32(const float complex *src, float *dst, uint32_t length);

void
dspm_sub_vf32_vf32_vf32(const float *src0,
                        const float *src1,
                        float *dst,
                        uint32_t length);

void
dspm_sub_vf32_vu32_vf32(const float *src0,
                        const uint32_t *src1,
                        float *dst,
                        uint32_t length);

void
dspm_clip_below_vf32_f32(float *srcdst,
                         float lb,
                         uint32_t length);

/* Note that this is not compensated summation */
float
dspm_sum_vf32(const float *src, uint32_t length);

/* Assumes d is non-zero */
void
dspm_div_vf32_f32(float *srcdst,
                  float d,
                  uint32_t length);

struct dspm_fft_init {
    /* The number of real values */
    uint32_t length;
    /* non-zero if inverse transform, zero if forward transform */
    uint32_t inverse;
    /*
    This structure can be subclassed as necessary for specific implementations
    (but then portability is compromised)
    */
};

/*
This is the "unpacked" real-only FFT implementation. That means the output
contains init->length/2 + 1 complex values where the first is purely real and
the last is purely imaginary.
*/
struct dspm_rfft_vf32_vz32_cfg *
dspm_rfft_vf32_vz32_cfg_new(struct dspm_fft_init *init);

void
dspm_rfft_vf32_vz32_cfg_free(struct dspm_rfft_vf32_vz32_cfg *cfg);

void
dspm_rfft_vf32_vz32(struct dspm_rfft_vf32_vz32_cfg *cfg,
                    float *time,
                    float complex *freq);

float
dspm_mean_vf32_f32(const float *src, uint32_t length);

void
dspm_rev_vu32(uint32_t *srcdst, uint32_t length);

void
dspm_floor_vf32_vu32(const float *src, uint32_t *dst, uint32_t length);

/* it is assumed that 0 <= i < length(src0) for all i in src1 and length(src1)
length(dst) = length */
void
dspm_lookup_vf32_vu32_vf32(
const float *src0,
const uint32_t *src1,
float *dst,
uint32_t length);

void
dspm_interp1d4p_vf32_vf32_vf32(const float *xi,
                               const float *y,
                               float *yi,
                               uint32_t N);

void
dspm_interp1d4p_vu24q8_vf32_vf32(const u24q8 *xi,
                                 const float *y,
                                 float *yi,
                                 uint32_t N);

void
dspm_interp1d4p_vu16q16_vf32_vf32(const u16q16 *xi,
                                  const float *y,
                                  float *yi,
                                  uint32_t N);

void
dspm_floor_vu24q8_vu32(const u24q8 *src,
                       uint32_t *dst,
                       uint32_t N);

void
dspm_floor_vu16q16_vu32(const u16q16 *src,
                       uint32_t *dst,
                       uint32_t N);

void
dspm_sub_vu24q8_vu32_vf32(const u24q8 *src0,
                          const uint32_t *src1,
                          float *dst,
                          uint32_t N);

/* computes dst[0] = initial_sum + src[0], dst[1] = initial_sum + src[0] +
src[1], etc. */
void
dspm_cumsum_vu16q16_u48q16_vu48q16(const u16q16 *src,
                                   u48q16 initial_sum,
                                   u48q16 *dst,
                                   uint32_t N);

dspm_cumsum_vs16q16_s48q16_vs48q16(const s16q16 *src,
                                   s48q16 initial_sum,
                                   s48q16 *dst,
                                   uint32_t N);

struct dspm_2dline_s48q16 {
    /* y = m*(x-x0) + b */
    s48q16 x0;
    s48q16 m;
    s48q16 b;
};

struct dspm_2dline_s48q16 dspm_2dline_s48q16_points(s48q16 x0,
                                                    s48q16 y0,
                                                    s48q16 x1,
                                                    s48q16 y1);

void
dspm_2dline_s48q16_lookup_vs48q16(const struct dspm_2dline_s48q16 *restrict line,
                                  s48q16 *restrict x,
                                  uint32_t N);

void
dspm_cvt_vu24q8_vf32(const u24q8 *src,
                     float *dst,
                     uint32_t N);

void
dspm_cvt_vu16q16_vf32(const u16q16 *src,
                     float *dst,
                     uint32_t N);

/* results are undefined if src is negative or greater than 2^16 */
void
dspm_cvt_vf32_vu16q16(const float *src,
                      u16q16 *dst,
                      uint32_t N);

/* For each element, puts the element or the constant, whichever is smaller. Can
be used to clip values to some upper boundary */
dspm_min_vu16q16_u16q16(u16q16 *srcdst, u16q16 c, uint32_t N);

/* For each element, puts the element or the constant, whichever is greater. Can
be used to clip values to some lower boundary */
dspm_max_vu16q16_u16q16(u16q16 *srcdst, u16q16 c, uint32_t N);

#endif /* DSP_MATH_H */
