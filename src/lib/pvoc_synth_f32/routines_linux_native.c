/*
Implementation of pvoc_synth that is not optimized for any particular
architecture but requires libraries generally available with every linux
distribution.
*/

#include "pvoc_synth_f32/routines_linux_native.h"
#include <stdlib.h>
#include <string.h>

#include "kiss_fftr.h"
#include "datastructures/ola_f32.h"

/* struct pvs_real_t's underlying implementation is an array of floats */

static struct pvs_real_t *
real_alloc(unsigned int length)
{
    return (struct pvs_real_t *)calloc(length,sizeof(float));
}

static void
real_free(struct pvs_real_t *s) { free(s); }

static void
real_memcpy(struct pvs_real_t *dst, const pvs_real_t *src, unsigned int N)
{
    memcpy(dst,src,sizeof(float)*N);
}

/* struct pvs_complex_t's underlying implementation is array of fftwf_complex */

/* The transform we will use is a transform of purely real values, so for a
transform of length N, we need only store N/2 + 1 complex values. */
static struct pvs_complex_t *
complex_alloc(unsigned int length)
{
    unsigned int N = length/2 + 1;
    fftwf_complex *ret = fftwf_malloc(sizeof(fftwf_complex) * N)
    if (!ret) { return ret; }
    memset(ret,0,sizeof(fftwf_complex) * N);
    return ret;
}

void
complex_free(struct pvs_complex_t *z) { fftwf_free((fftwf_complex*)z); }

/*
in the underlying implementation of struct pvs_dft_t we store the "plans" for
the forward and backward fourier transforms
*/
struct pvs_dft_t {
    kiss_fftr_cfg forward_dft;
    kiss_fftr_cfg inverse_dft;
};

static void
dft_free(struct pvs_dft_t *dft)
{
    if (!dft) { return; }
    if (dft->forward_dft) { kiss_fft_free(dft->forward_dft); }
    if (dft->inverse_dft) { kiss_fft_free(dft->inverse_dft); }
    free(dft);
}

static struct pvs_dft_t *
dft_alloc(pvs_dft_init_t *init)
{
    if (!init) { return NULL; }
    struct pvs_dft_t *ret = calloc(1,sizeof(struct pvs_dft_));
    if (!ret) { goto fail; }
    ret->forward_dft = kiss_fftr_alloc(init->N,0,NULL,NULL);
    if (!ret->forward_dft) { goto fail; }
    ret->inverse_dft = kiss_fftr_alloc(init->N,1,NULL,NULL);
    if (!ret->inverse_dft) { goto fail; }
    return ret
fail:
    dft_free(ret);
    return NULL;
}

static void
dft_window_scale(struct pvs_dft_window_scale_t *dws)
{
    /* Simply scale synthesis window by 1/N */
    float scalar = 1./dws->window_length;
    float *window = (float*)dws->synthesis_window;
    unsigned int n;
    for (n = 0; n < dws->window_length; n++) {
        window[n] *= scalar;
    }
}

/* perform forward DFT */
static void dft_forward (struct pvs_dft_t * dft, const struct pvs_real_t * a,
                     struct pvs_complex_t * b)
{
    kiss_fftr(dft->forward_dft,(const kiss_fft_scalar*)a,(kiss_fft_cpx*)b);
}

/* perform inverse DFT */
static void dft_inverse (struct pvs_dft_t * dft,
                     const struct pvs_complex_t * a,
                     struct pvs_real_t * b)
{
    kiss_fftri(dft->inverse_dft,(const kiss_fft_cpx*)a,(kiss_fft_scalar*)b);
}

#define _ARRAY_MULTIPLY(a,b,t1,t2,l)\
    t1 a_ = (t1)a;\
    t2 b_ = (t2)b;\
    while (l--) { *a_++ *= *b_++; }

#define _ARRAY_DIVIDE(a,b,t1,t2,l)\
    t1 a_ = (t1)a;\
    t2 b_ = (t2)b;\
    while (l--) { *a_++ /= *b_++; }

/* multiply 2 complex arrays, result goes in first array (i.e., a *= b) */
static void complex_complex_mult (struct pvs_complex_t * a,
                              const struct pvs_complex_t * b,
                              unsigned int length)
{
    _ARRAY_MULTIPLY(a,b,complex float*,const complex float*,length);
}

/* multiply complex and real arrays, result goes in first array (i.e., a *= b) */
static void complex_real_mult (struct pvs_complex_t * a,
                           const struct pvs_real_t * b,
                           unsigned int length)
{
    _ARRAY_MULTIPLY(a,b,complex float*,const float*,length);
}

/* multiply 2 real arrays, result goes in first array (i.e., a *= b) */
static void real_real_mult (struct pvs_real_t * a,
                        const struct pvs_real_t * b, unsigned int length)
{
    _ARRAY_MULTIPLY(a,b,float*,const float*,length);
}

/* multiply 2 real arrays, result goes in new array (i.e., c = a * b) */
static void real_real_cpymult (const struct pvs_real_t * a,
                        const struct pvs_real_t * b,
                        struct pvs_real_t * c,
                        unsigned int length)
{
    const float *a_ = (const float*)a,
                 b_ = (const float *)b;
    float *c_ = (float*)c;
    while (length--) {
        *c_++ = *a_++ * *b_++;
    }
}

/* divide 2 complex arrays, result goes in first array (i.e., a /= b) */
static void complex_complex_div (pvs_complex_t * a, const pvs_complex_t * b,
                             unsigned int length)
{
    _ARRAY_DIVIDE(a,b,complex float *, const complex float *,length);
}

/* put absolute value (modulus) of the complex values in an array. */
static void complex_abs (const pvs_complex_t * src, pvs_real_t * dst,
                     unsigned int length)
{
    const complex float *src_ = (const complex float*)src;
    float *dst_ = (float *)dst;
    while (length--) {
        *dst_++ = cabs(*src++);
    }
}

/* a table of functions implementing the required functionality */
static struct pvs_func_table_t func_table =
{
  .dstructs = {
    .real_alloc = real_alloc,
    .real_free = real_free,
    .real_memcpy = real_memcpy,
    .complex_alloc = complex_alloc,
    .complex_free = complex_free,
    .dft_free = dft_free,
    .dft_alloc = dft_alloc,
    .dft_window_scale = dft_window_scale,
    .ola_alloc = (struct pvs_ola_t *(*) (pvs_ola_init_t *))ola_f32_new,
    .ola_free = (void (*) (struct pvs_ola_t *))ola_f32_free,
    .ola_sum_in = (void (*) (struct pvs_ola_t *, const struct pvs_real_t *))ola_f32_sum_in,
    .ola_shift_out = (const struct pvs_real_t *(*) (struct pvs_ola_t *))ola_f32_shift_out,
    .dft_alloc = dft_alloc,
    .dft_free = dft_free
  },
  .math = {
    .complex_complex_mult = complex_complex_mult,
    .complex_real_mult = complex_real_mult,
    .real_real_mult = real_real_mult,
    .real_real_cpymult = real_real_cpymult,
    .complex_complex_div = complex_complex_div,
    .complex_abs = complex_abs,
    .dft_forward = dft_forward,
    .dft_inverse = dft_inverse
  } math;
};


struct pvs_init_t
pvs_f32_init_new (void)
{
    static struct pvs_init_t ret = {
        .func_table = func_table;
    };
    return ret;
}
