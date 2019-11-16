#ifndef LATTICE_FILTER_F32_H
#define LATTICE_FILTER_F32_H 

struct lattice_filter_f32;

struct lattice_filter_f32_init {
    unsigned int P;
};

enum lattice_filter_f32_opts {
    lattice_filter_f32_NONE = 0,
    /*
    The P reflection coefficients array R and the zeros converted to C
    coefficients are in fact an array of size P*N with R[0] to R[P-1] the
    coefficients used for the first sample, R[P] to R[2*P-1] the coefficients
    used for the second sample, etc. (the same goes for the C coefficients
    except that there are P+1 c coefficients per time step).
    */
    lattice_filter_f32_VARYING_R_C = (1 << 0),
};

struct lattice_filter_f32_proc {
    const float *in;
    float *out;
    unsigned int N;
    /* P reflection coefficients. */
    const float *R;
    /* P+1 "C" coefficients. */
    const float *C;
    const enum lattice_filter_f32_opts opts;
};

void lattice_filter_f32_free(struct lattice_filter_f32 *l);
struct lattice_filter_f32 * lattice_filter_f32_new( struct lattice_filter_f32_init *init);
void lattice_filter_f32_process( struct lattice_filter_f32 *l, struct lattice_filter_f32_proc *p);

#endif /* LATTICE_FILTER_F32_H */
