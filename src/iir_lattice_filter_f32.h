#ifndef IIR_LATTICE_FILTER_F32_H
#define IIR_LATTICE_FILTER_F32_H 

struct iir_lattice_filter_f32;

struct iir_lattice_filter_f32_init {
    unsigned int P;
};

struct iir_lattice_filter_f32_proc {
    const float *in;
    float *out;
    unsigned int N;
    /* P reflection coefficients. */
    const float *R;
    /* The first and only feedforward coefficient */
    const float b0;
};

void iir_lattice_filter_f32_free(struct iir_lattice_filter_f32 *l);
struct iir_lattice_filter_f32 * iir_lattice_filter_f32_new( struct iir_lattice_filter_f32_init *init);
void iir_lattice_filter_f32_process( struct iir_lattice_filter_f32 *l, struct iir_lattice_filter_f32_proc *p);

#endif /* IIR_LATTICE_FILTER_F32_H */
