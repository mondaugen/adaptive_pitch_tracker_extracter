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
    /* Allocate array of reals */
    struct pvs_real_t *(*real_alloc) (unsigned int length);
    /* free array of reals */
    void (*real_free) (struct pvs_real_t *);
    /* Copy array of reals */
    void (*real_memcpy)(struct pvs_real_t *, const struct pvs_real_t *, unsigned int);
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
    /* Allocate and initialize DFT auxiliary structure */
    struct pvs_dft_t *(*dft_alloc) (struct pvs_dft_init_t *);
    /* Free DFT auxiliary structure */
    void (*dft_free) (struct pvs_dft_t *);
    /* A function that scales the analysis and / or synthesis windows so that
    the forward followed by inverse transform of the windowed signal is not
    scaled by any constant (for example, some implementations will omit the 1/N
    on the inverse transform so that the composition of the forward and inverse
    transforms produces the original sequence but divided by N) */
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

struct pvs_init_t
{
  /* The resulting pvs_t makes copies of the analysis and synthesis windows */ 
  /* The signal is localized in time by multiplying this window (assumed 0
     outside of this array of length window_length). */
  const pvs_real_t *analysis_window;
  /* The output is multiplied by this synthesis window */
  const pvs_real_t *synthesis_window;
  /* Length of analysis and synthesis windows */
  unsigned int window_length;
  /* The number of samples between the beginnings of the 2 analysis windows */
  unsigned int hop_size;
  /* A function that is passed the auxilary data structure, and a
     pvs_f32_sample_lookup_t structure, which then fills this with the number of
     samples */
  struct pvs_real_t *(*get_samples) (void *aux,
                                     struct pvs_real_sample_lookup_t * info);
  /* Auxiliary structure for get_samples */
  void *get_samples_aux;
  /* Table of functions implementing functionality */
  struct pvs_func_table_t * const func_table;
};

struct pvs_t;

void pvs_process(pvs_t *pvs, struct pvs_real_t *output, int input_time);
void pvs_free(struct pvs_t *pvs);
struct pvs_t * pvs_new(pvs_init_t *init);

#endif /* PVOC_SYNTH_H */
