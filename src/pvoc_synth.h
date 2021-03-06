#ifndef PVOC_SYNTH_H
#define PVOC_SYNTH_H

/*
An implementation must subclass these types and provide function
implementations by passing a vector table.
*/

/* array of reals */
struct pvs_real_t;

/* Functions implementing get_samples use this structure */
struct pvs_real_sample_lookup_t {
    /* the get_samples function should return samples for any value of first_sample_index */
    int first_sample_index;
    unsigned int n_samples;
    const struct pvs_real_t *samples;
};

/* array of complex numbers */
struct pvs_complex_t;

/* DFT initialization structure.  */
struct pvs_dft_init_t
{
  /* the size of the transform */
  unsigned int N;
};

/* DFT auxiliary stuff */
struct pvs_dft_t;

/* Window scale function argument */
struct pvs_dft_window_scale_t {
  struct pvs_dft_t *dft;
  struct pvs_real_t *analysis_window;
  struct pvs_real_t *synthesis_window;
  /* Length of analysis and synthesis windows */
  unsigned int window_length;
};

/* Overlap-and-add buffer initialization structure. The structure is designed to
work with fixed summing-in and shifting-out lengths */
struct pvs_ola_init_t
{
  /* This is usually equal to the window length of the transform */
  unsigned int sum_in_length;
  /* This is usually equal to the hop size or equivalently the audio
     processing block size */
  unsigned int shift_out_length;
};

/* Overlap-and-add buffer for reals */
struct pvs_ola_t;

/* a table of functions implementing the required functionality */
struct pvs_func_table_t
{
  struct
  {
    /* NOTE: The allocation functions must return zeroed memory (e.g., by using calloc)*/
    /* Allocate array of reals */
    struct pvs_real_t *(*real_alloc) (unsigned int length);
    /* free array of reals */
    void (*real_free) (struct pvs_real_t *);
    /* Copy array of reals */
    void (*real_memcpy)(struct pvs_real_t *, const struct pvs_real_t *, unsigned int);
    /*
    Get a pointer to an array but starting at offset (similar to x[offset:]
    in Python)
    */
    struct pvs_real_t *(*real_offset) (struct pvs_real_t *r, unsigned int offset);
    /* Allocate array of complex */
    struct pvs_complex_t *(*complex_alloc) (unsigned int length);
    /* free array of complex */
    void (*complex_free) (struct pvs_complex_t *);
    /* Allocate real overlap-and-add buffer */
    struct pvs_ola_t *(*ola_alloc) (struct pvs_ola_init_t *);
    /* Free real overlap-and-add buffer */
    void (*ola_free) (struct pvs_ola_t *);
    /* Sum in sum_in_length values into a pvs_ola_t and shift out
    shift_out_length values from a pvs_ola_t into a pvs_real_t.  This returns a
    pointer to contiguous memory holding hop_size real values and advances the
    pvs_ola_t's internal pointer by hop_size.  */
    struct pvs_real_t * (*ola_sum_in_and_shift_out) (struct pvs_ola_t *, const struct pvs_real_t *);
    /* Allocate and initialize DFT auxiliary structure */
    struct pvs_dft_t *(*dft_alloc) (struct pvs_dft_init_t *);
    /* Free DFT auxiliary structure */
    void (*dft_free) (struct pvs_dft_t *);
    /* A function that scales the analysis and / or synthesis windows so that
    the forward followed by inverse transform of the windowed signal is not
    scaled by any constant (for example, some implementations will omit the 1/N
    on the inverse transform so that the composition of the forward and inverse
    transforms produces the original sequence but multiplied by N) */
    void (*dft_window_scale) (struct pvs_dft_window_scale_t *);
  } dstructs;
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
    /* multiply 2 real arrays, summing into other array */
    void (*real_real_add_product) (const struct pvs_real_t * a,
                            const struct pvs_real_t * b,
                            struct pvs_real_t * c,
                            unsigned int length);
    /* replace an array with its reciprocal (i.e. a = 1/a) */
    void (*real_reciprocal) (struct pvs_real_t * a, unsigned int length);
    /* returns non-zero if the array contains 0, 0 otherwise */
    int (*real_contains_zero) (const struct pvs_real_t *a, unsigned int length);
    /* divide 2 complex arrays, result goes in first array (i.e., a /= b) */
    void (*complex_complex_div) (struct pvs_complex_t * a,
                                 const struct pvs_complex_t * b,
                                 unsigned int length);
    /* put absolute value (modulus) of the complex values in an array. */
    void (*complex_abs) (const struct pvs_complex_t * src, struct pvs_real_t * dst,
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

struct pvs_init_t
{
    struct pvs_user_init_t {
        /* The resulting pvs_t makes copies of the analysis and synthesis windows */ 
        /* The signal is localized in time by multiplying this window (assumed 0
           outside of this array of length window_length). */
        const struct pvs_real_t *analysis_window;
        /* The output is multiplied by this synthesis window */
        const struct pvs_real_t *synthesis_window;
        /* Length of analysis and synthesis windows */
        unsigned int window_length;
        /* The number of samples between the beginnings of the 2 analysis windows */
        unsigned int hop_size;
        /* A function that is passed the auxilary data structure, and a
           pvs_f32_sample_lookup_t structure, which then fills this with the number of
           samples. 
           NOTE: For some pathological signals (which are most likely not
           encountered in practice) you can get a power spectrum where some bins
           have exactly 0 magnitude. This is problematic for this algorithm,
           which finds the quotient of 2 complex numbers. In practice this
           shouldn't happen if a small amount of dither is added to the signal,
           but it could happen. In that case, somewhere in the system a check
           should be made for NaN, inf and the like. */
        void (*get_samples) (void *aux,
                struct pvs_real_sample_lookup_t * info);
        /* Auxiliary structure for get_samples */
        void *get_samples_aux;
    } user;
  /* Table of functions implementing functionality */
  struct pvs_func_table_t * const func_table;
};

struct pvs_t;

const struct pvs_real_t * pvs_process(struct pvs_t *pvs, int input_time);
void pvs_free(struct pvs_t *pvs);
struct pvs_t * pvs_new(struct pvs_init_t *init);

#endif /* PVOC_SYNTH_H */
