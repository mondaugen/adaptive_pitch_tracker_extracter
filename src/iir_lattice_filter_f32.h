#ifndef IIR_LATTICE_FILTER_F32_H
#define IIR_LATTICE_FILTER_F32_H 

struct iir_lattice_filter_f32;

struct iir_lattice_filter_f32_init {
    unsigned int P;
};

enum iir_lattice_filter_f32_opts {
    iir_lattice_filter_f32_NONE = 0,
    /*
    The P reflection coefficients array R is in fact an array of size P*N
    with R[0] to R[P-1] the coefficients used for the first sample, R[P] to
    R[2*P-1] the coefficients used for the second sample, etc.
    */
    iir_lattice_filter_f32_VARYING_R = (1 << 0),
};

struct iir_lattice_filter_f32_proc {
    const float *in;
    float *out;
    unsigned int N;
    /* P reflection coefficients. */
    const float *R;
    /* The first and only feedforward coefficient */
    const float b0;
    const enum iir_lattice_filter_f32_opts opts;
};

void iir_lattice_filter_f32_free(struct iir_lattice_filter_f32 *l);
struct iir_lattice_filter_f32 * iir_lattice_filter_f32_new( struct iir_lattice_filter_f32_init *init);
void iir_lattice_filter_f32_process( struct iir_lattice_filter_f32 *l, struct iir_lattice_filter_f32_proc *p);

#endif /* IIR_LATTICE_FILTER_F32_H */
