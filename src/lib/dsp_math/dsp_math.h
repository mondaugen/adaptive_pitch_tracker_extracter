#ifndef DSP_MATH_H
#define DSP_MATH_H 

#include <complex.h>

/*
DSP math routines
The naming convention here is similar to Intel's IPP.
functions begin with dspm_
then a short description of the operation (e.g., dspm_add)
then the types, prepended optionally with v if a vector is accepted as an
argument (a pointer), followed by the type name (e.g., f32 for 32-bit float).
so the function dspm_add_vf32_vf32(float *srcdst, const float* src, unsigned int
length) computes srcdst[n] += src[n] for n in [0,length) and the function
dspm_add_vf32_vf32_vf32(const float *src0, const float *src1, float *dst,
unsigned int length) computs dst[n] = src0[n] + src1[n] for n in [0,length).

For operations requiring vectors, calling these on vectors of length zero or
NULL pointer arguments, the results are undefined.

*/

/* Fast floor(log2(x)) */
static inline int
dspm_fast_floor_log2_f32(float x)
{
    unsigned int *ux = (unsigned int*)&x;
    /* Extracts exponent from IEEE 754 floating point number */
    int exp = ((*ux>>23)&0xff) - 127;
    return exp;
}

/* Fast approximate fractional part of log2(x) */
static inline float
dspm_fast_log2_aprox_frac_f32(float x)
{
    unsigned int *ux = (unsigned int*)&x;
    /* (*ux&0x7fffff)**(2**-23) */
    float frac = (*ux&0x7fffff)*1.1920928955078125e-07;
    return frac;
}
    
void
dspm_mul_vf32_vf32_vf32(const float *src0,
                        const float *src1,
                        float *dst,
                        unsigned int length);

void
dspm_mul_vf32_f32(float *srcdst,
                  float src1,
                  unsigned int length);

void
dspm_neg_vf32(float *srcdst,
              unsigned int length);

float
dspm_max_vf32(const float *src,
              unsigned int length);

void
dspm_abs_vf32(float *srcdst, unsigned int length);

void
dspm_abs_vz32(float complex *srcdst, unsigned int length);

void
dspm_abs_vz32_vf32(const float complex *src, float *dst, unsigned int length);

void
dspm_sub_vf32_vf32_vf32(const float *src0,
                        const float *src1,
                        float *dst,
                        unsigned int length);

void
dspm_clip_below_vf32_f32(float *srcdst,
                         float lb,
                         unsigned int length);

/* Note that this is not compensated summation */
float
dspm_sum_vf32(const float *src, unsigned int length);

/* Assumes d is non-zero */
void
dspm_div_vf32_f32(float *srcdst
                  float d,
                  unsigned int length);

struct dspm_fft_init {
    /* The number of real values */
    unsigned int length;
    /* non-zero if inverse transform, zero if forward transform */
    unsigned int inverse;
    /*
    This structure can be subclassed as necessary for specific implementations
    (but then portability is compromised).
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
                    const float *src,
                    float complex *dst);

#endif /* DSP_MATH_H */
