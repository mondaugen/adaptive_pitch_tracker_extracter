#include <stdlib.h>
#include "gal_f32.h"

#include <stdio.h>

#define conj(x) (x)

struct gal_f32 {
    unsigned int P; /* Filter order */
    /* The backward errors at each stage for this time step, length P + 1 */
    float *Eb_n_j;
    /* The backward errors at each stage for the previous time step, length P + 1 */
    float *Eb_n_1_j;
    /* The previous reflection coefficients, length P */
    float *R;
    /* The previous normalization coefficients */
    float *D;
};

void
gal_f32_free(struct gal_f32 *g)
{
    if (g->Eb_n_j) { free(g->Eb_n_j); }
    if (g->Eb_n_1_j) { free(g->Eb_n_1_j); }
    if (g->R) { free(g->R); }
    if (g->D) { free(g->D); }
    free(g);
}

struct gal_f32 *
gal_f32_new(struct gal_f32_init *init)
{
    if (init->P < 1) { return NULL; }
    struct gal_f32 *ret = calloc(1,sizeof(struct gal_f32));
    if (!ret) { return NULL; }
    ret->P = init->P;
    ret->Eb_n_j = calloc(ret->P + 1,sizeof(float));
    if (!ret->Eb_n_j) { goto fail; }
    ret->Eb_n_1_j = calloc(ret->P + 1, sizeof(float));
    if (!ret->Eb_n_1_j) { goto fail; }
    ret->R = calloc(ret->P, sizeof(float));
    if (!ret->R) { goto fail; }
    ret->D = calloc(ret->P, sizeof(float));
    if (!ret->D) { goto fail; }
    /* Initialize to 1s */
    int p;
    for (p = 0 ; p < ret->P; p++) {
        ret->D[p] = 1;
    }
    return ret;
fail:
    gal_f32_free(ret);
    return NULL;
}

void
gal_f32_process(struct gal_f32 *g, struct gal_f32_proc *p)
{
    unsigned int n, j;
    /* The forward errors for the past-stage and current stage at time step n */
    float Ef_j_1_n, Ef_j_n, *tmp;
    for (n = 0; n < p->N; n++) {
        Ef_j_1_n = p->x_in[n];
        g->Eb_n_j[0] = p->x_in[n];
        for (j = 1; j <= g->P; j++) {
            Ef_j_n = Ef_j_1_n + g->R[j-1] * g->Eb_n_1_j[j-1];
            g->Eb_n_j[j] = g->Eb_n_1_j[j-1] + conj(g->R[j-1]) * Ef_j_1_n;
            if (p->opt & gal_f32_opts_ESTIMATE_MU) {
                p->D[g->P*n + j - 1] = g->D[j-1] = (1 - p->alpha) * g->D[j-1] +
                    p->alpha * (Ef_j_n * Ef_j_n + g->Eb_n_j[j] * g->Eb_n_j[j]);
            }
            /* otherwise D remains its old value */
            /* Update and store R
               TODO: use pointer arithmetic to save address computation */
            g->R[j-1] = p->R[n*g->P + j - 1] = g->R[j-1] - p->beta * (Ef_j_n *
                conj(g->Eb_n_1_j[j-1]) + conj(g->Eb_n_j[j]) * Ef_j_1_n) / g->D[j-1];
            Ef_j_1_n = Ef_j_n;
        }
        /* Store final stage */
        p->Ef[n] = Ef_j_n;
        p->Eb[n] = g->Eb_n_j[g->P - 1];
        /* Swap to make current backward errors the previous ones */
        tmp = g->Eb_n_j;
        g->Eb_n_j = g->Eb_n_1_j;
        g->Eb_n_1_j = tmp;
    }
}
