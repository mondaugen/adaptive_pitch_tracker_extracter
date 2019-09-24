#ifndef PVOC_SYNTH_F32_H
#define PVOC_SYNTH_F32_H 

/*
A struct holding stuff for doing a DFT, to be defined in an implementation.
*/
struct pvs_f32_dft_t;

/*
Its initialization structure. This can be subclassed as needed if
implementations need more information.
*/
struct pvs_f32_dft_init_t {
    /* the size of the transform */
    unsigned int N;
};

/*
A struct holding 64 bit complex values (2 packed f32s)
*/
struct pvs_z32_array_t {
    /* The number of complex values in the array */
    unsigned int length;
};

/* Function allocating an array of complex64 values. */
struct pvs_z32_array_new(unsigned int length);

/* Function freeing array of complex64 values. */
void
pvs_z32_array_free(struct pvs_z32_array_t *a);

/*
Function multiplying 2 z32 arrays, result goes in first array (i.e., a *= b)
*/
void
pvs_z32_z32_mult(pvs_z32_array_t *a, const pvs_z32_array_t *b);

/*
Function dividing 2 z32 arrays, result goes in first array (i.e., a /= b)
*/
void
pvs_z32_z32_div(pvs_z32_array_t *a, const pvs_z32_array_t *b);

/*
Function putting the absolute value (modulus) of the complex values in an array.
dst must have length equal to src->length.
*/
void
pvs_z32_abs(const pvs_z32_array_t *src, float *dst);

/*
Function multiplying z32 with f32 array, result goes in first array
(i.e., a *= b). b must have length equal to a->length.
*/
void
pvs_z32_f32_mult(pvs_z32_array_t *a, const float *b);

/*
Performs forward DFT on a, putting the result in b. a must have length equal to b->length.
*/
void
pvs_f32_dft_forward(struct pvs_f32_dft_t *dft, const float *a, pvs_z32_array_t *b);

/*
Performs inverse DFT on a, putting result in b. b must have length equal to a->length.
*/
void
pvs_f32_dft_inverse(struct pvs_f32_dft_t *dft, const pvs_z32_array_t *a, float *b);

/* Functions implementing get_samples use this structure */
struct pvs_f32_sample_lookup_t {
    unsigned int first_sample_index;
    unsigned int n_samples;
    float *samples;
};

#endif /* PVOC_SYNTH_F32_H */
