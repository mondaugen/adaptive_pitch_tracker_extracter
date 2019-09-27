/*
Implementation of pvoc_synth that is not optimized for any particular
architecture but requires libraries generally available with every linux
distribution.
*/

#include "kiss_fftr.h"
#include "pvoc_synth.h"
#include "datastructures/ola_f32.h"

/* struct pvs_real_t's underlying implementation is an array of floats */

static struct pvs_real_t *
real_alloc(unsigned int length)
{
    return (struct pvs_real_t *)calloc(length,sizeof(float));
}

static void
real_free(struct pvs_real_t *s) { free(s); }

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

/* a table of functions implementing the required functionality */
static struct pvs_func_table_t func_table =
{
  .dstructs = {
    /* Allocate array of reals */
    struct pvs_real_t *(*real_alloc) (unsigned int length);
    /* free array of reals */
    void (*real_free) (struct pvs_real_t *);
    /* Allocate array of complex */
    struct pvs_complex_t (*complex_alloc) (unsigned int length);
    /* free array of complex */
    void (*complex_free) (struct pvs_complex_t *);
    /* Allocate real overlap-and-add buffer */
    struct pvs_ola_t *(*ola_alloc) (pvs_ola_init_t *);
    /* Free real overlap-and-add buffer */
    void (*ola_free) (struct pvs_ola_t *);
    /* Sum in sum_in_length values into a pvs_ola_t */
    void (*ola_sum_in) (struct pvs_ola_t *, const struct pvs_real_t *);
    /*
    Shift out shift_out_length values from a pvs_ola_t into a pvs_real_t.
    This returns a pointer to contiguous memory holding hop_size real values and
    advances the pvs_ola_t's internal pointer by hop_size.
    */
    const struct pvs_real_t *(*ola_shift_out) (struct pvs_ola_t *);
    /* Allocate DFT auxiliary structure */
    struct pvs_dft_t *(*dft_alloc) (pvs_dft_init_t *);
    /* Free DFT auxiliary structure */
    void (*dft_free) (pvs_dft_t *);
  },
  struct
  {
    /* multiply 2 complex arrays, result goes in first array (i.e., a *= b) */
    void (*complex_complex_mult) (struct pvs_complex_t * a,
                                  const struct pvs_complex_t * b,
                                  unsigned int length);
    /* multiply complex and real arrays, result goes in first array (i.e., a *= b) */
    void (*complex_real_mult) (struct pvs_complex_t * a,
                               const struct pvs_real_t * b,
                               unsigned int length);
    /* multiply 2 real arrays, result goes in first array (i.e., a *= b) */
    void (*real_real_mult) (struct pvs_real_t * a,
                            const struct pvs_real_t * b, unsigned int length);
    /* multiply 2 real arrays, result goes in new array (i.e., c = a * b) */
    void (*real_real_cpymult) (const struct pvs_real_t * a,
                            const struct pvs_real_t * b,
                            struct pvs_real_t * c,
                            unsigned int length);
    /* divide 2 complex arrays, result goes in first array (i.e., a /= b) */
    void (*complex_complex_div) (pvs_complex_t * a, const pvs_complex_t * b,
                                 unsigned int length);
    /* put absolute value (modulus) of the complex values in an array. */
    void (*complex_abs) (const pvs_complex_t * src, pvs_real_t * dst,
                         unsigned int length);
    /* perform forward DFT */
    void (*dft_forward) (struct pvs_dft_t * dft, const struct pvs_real_t * a,
                         struct pvs_complex_t * b);
    /* perform inverse DFT */
    void (*dft_inverse) (struct pvs_dft_t * dft,
                         const struct pvs_complex_t * a,
                         struct pvs_real_t * b);
  } math;
};

