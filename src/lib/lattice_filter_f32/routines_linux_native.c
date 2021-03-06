/*

Lattice Filter

The notation used here is consistent with the notation in "Statistical Signal
Processing and Modeling" by Monson Hayes. See chapter 6.

For users:
The structure implemented has input values entering on the left side of the
structure and output values exiting on the right side of the structure.
The reflection coefficients are stored so the first reflection coefficient is
the one on the right-most side of the structure. Furthermore, referring to the
book above, R[0] corresponds to Gamma_1 (because Gamma_0 is never used or
referenced).
Basically, if you store the feedback coefficients so that the first one is the
y[n] coefficient, the second the y[n-1] coefficient, etc., and you use the
step-down recursion as described in the book to obtain the reflection
coefficients, these resulting reflection coefficients can be used in the same
order for this filter implementation.

The C coefficients can be computed as described on p. 306 in the Hayes book.

*/

#include <stdlib.h>
#include "lattice_filter_f32.h"

struct lattice_filter_f32 {
    /* Order and number of reflection coefficients */
    unsigned int P;
    /* The P e- values i.e., the backward errors */
    float *e_m;
    /* The P e+ values i.e., the forward errors */
    float *e_p;
    /* The P past e+ values */
    float *e_m_z1;
};

void
lattice_filter_f32_free(struct lattice_filter_f32 *l)
{
    if (!l) { return; }
    if (l->e_m) { free(l->e_m); }
    if (l->e_p) { free(l->e_p); }
    if (l->e_m_z1) { free(l->e_m_z1); }
    free(l);
}

struct lattice_filter_f32 *
lattice_filter_f32_new(
    struct lattice_filter_f32_init *init)
{
    if (!init) { return NULL; }
    struct lattice_filter_f32 *ret = calloc(1,sizeof(struct lattice_filter_f32));
    if (!ret) { return NULL; }
    ret->P = init->P;
    ret->e_m = calloc(ret->P+1,sizeof(float));
    if (!ret->e_m) { goto fail; }
    ret->e_p = calloc(ret->P+1,sizeof(float));
    if (!ret->e_p) { goto fail; }
    ret->e_m_z1 = calloc(ret->P+1,sizeof(float));
    if (!ret->e_m_z1) { goto fail; }
    return ret;
fail:
    lattice_filter_f32_free(ret);
    return NULL;
}

void
lattice_filter_f32_process(
    struct lattice_filter_f32 *l,
    struct lattice_filter_f32_proc *p)
{
    float *out = p->out, *e_m, *e_p, *e_m_z1, *tmp;
    const float *in = p->in, *R, *C, *p_R = p->R, *p_C = p->C;
    unsigned int N = p->N, P, p_R_inc = 0, p_C_inc = 0;
    if (p->opts & lattice_filter_f32_VARYING_R_C) {
        p_R_inc = l->P;
        p_C_inc = l->P+1;
    }
    while (N--) {
        P = l->P;
        e_p = l->e_p + P;
        *e_p = *in;
        e_p--;
        R = p_R + P - 1;
        e_m_z1 = l->e_m_z1 + P - 1;
        while (P--) {
            *e_p = *(e_p+1) - *R * *e_m_z1;
            e_p--;
            R--;
            e_m_z1--;
        }
        /* Propagate the backward errors */
        e_m = l->e_m;
        e_m_z1 = l->e_m_z1;
        e_p = l->e_p;
        e_m = l->e_m;
        *e_m = *e_p;
        R = p_R;
        P = l->P;
        while (P--) {
            *(e_m + 1) = *e_m_z1 + *R * *e_p;
            e_m++;
            e_m_z1++;
            e_p++;
            R++;
        }
        /* Add up backwards errors weighted by C coefficients */
        C = p_C;
        e_m = l->e_m;
        *out = *e_m++ * *C++;
        P = l->P;
        while (P--) {
            *out += *e_m++ * *C++;
        }
        /* Swap in order to delay */
        tmp = l->e_m;
        l->e_m = l->e_m_z1;
        l->e_m_z1 = tmp;
        in++;
        out++;
        p_R += p_R_inc;
        p_C += p_C_inc;
    }
}
