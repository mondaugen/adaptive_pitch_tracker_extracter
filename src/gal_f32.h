#ifndef GAL_F32_H
#define GAL_F32_H

struct gal_f32;

struct gal_f32_init {
    unsigned int P;
};

enum gal_f32_opts {
    gal_f32_opts_NONE = 0,
    /* If the ESTIMATE_MU option is passed then the filter will try and estimate
    the mu coefficients using the beta and lambda provided (see p. 528 Hayes). */
    gal_f32_opts_ESTIMATE_MU = (1 << 0),
};

struct gal_f32_proc {
    /* The signal to be modelled, length N */
    const float *x_in;
    /* The forward error signal, i.e., the output of the adaptive FIR filter,
    length N */
    float *ep;
    /* The backward error signal, length N */
    float *em;
    /* The estimated reflection coefficients, length N*P */
    float *R;
    /* The gradient descent coefficients, length P. If opt has ESTIMATE_MU set,
    then memory for these must still be provided but they will be estimated from
    beta and lambda below. Otherwise beta and lambda are ignored */
    float *mu;
    float beta;
    float lambda;
    enum gal_f32_opts opt;
    unsigned int N;
};

extern void
gal_f32_free(struct gal_f32 *g);

extern struct gal_f32 *
gal_f32_new(struct gal_f32_init *init);

extern void
gal_f32_process(struct gal_f32 *g, struct gal_f32_proc *p);

#endif /* GAL_F32_H */
