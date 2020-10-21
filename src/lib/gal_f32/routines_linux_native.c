#include <stdlib.h>
#include "gal_f32.h"

#include <stdio.h>

#define conj(x) (x)

struct gal_f32 {
    unsigned int P; /* Filter order */
    /* The backward errors at each stage for this time step, length P + 1 */
    float *em_n_j;
    /* The backward errors at each stage for the previous time step, length P + 1 */
    float *em_n_1_j;
    /* The current reflection coefficients, length P */
    float *R;
};

void
gal_f32_free(struct gal_f32 *g)
{
    if (g->em_n_j) { free(g->em_n_j); }
    if (g->em_n_1_j) { free(g->em_n_1_j); }
    if (g->R) { free(g->R); }
    free(g);
}

struct gal_f32 *
gal_f32_new(struct gal_f32_init *init)
{
    if (init->P < 1) { return NULL; }
    struct gal_f32 *ret = calloc(1,sizeof(struct gal_f32));
    if (!ret) { return NULL; }
    ret->P = init->P;
    ret->em_n_j = calloc(ret->P + 1,sizeof(float));
    if (!ret->em_n_j) { goto fail; }
    ret->em_n_1_j = calloc(ret->P + 1, sizeof(float));
    if (!ret->em_n_1_j) { goto fail; }
    ret->R = calloc(ret->P, sizeof(float));
    if (!ret->R) { goto fail; }
    return ret;
fail:
    gal_f32_free(ret);
    return NULL;
}

void
gal_f32_process(struct gal_f32 *g, struct gal_f32_proc *p)
{
    unsigned int n, j;
    float ep_j_1_n, ep_j_n, *tmp;
    for (n = 0; n < p->N; n++) {
        ep_j_1_n = p->x_in[n];
        g->em_n_j[0] = p->x_in[n];
        for (j = 1; j <= g->P; j++) {
            ep_j_n = ep_j_1_n + g->R[j-1] * g->em_n_1_j[j-1];
            g->em_n_j[j] = g->em_n_1_j[j-1] + conj(g->R[j-1]) * ep_j_1_n;
            /* Update and store R
               TODO: use pointer arithmetic to save address computation */
            g->R[j-1] = p->R[n*g->P + j - 1] = g->R[j-1] - p->mu[j-1] * (ep_j_n *
                conj(g->em_n_1_j[j]) + conj(g->em_n_j[j]) * ep_j_1_n);
            if (p->opt |= gal_f32_opts_ESTIMATE_MU) {
                p->D[j-1] = p->lambda * p->D[j-1] + (1 - p->lambda) * (
                    ep_j_n * ep_j_n + g->em_n_1_j[j] * g->em_n_1_j[j]
                );
                /* Avoid divide by 0 */
                p->mu[j - 1] = p->beta / (p->D[j-1] + 1e-6);
                //printf("%f ", p->D[j-1]);
            }
            ep_j_1_n = ep_j_n;
        }
        /* Store final stage */
        p->ep[n] = ep_j_n;
        p->em[n] = g->em_n_j[g->P - 1];
        /* Swap to make current backward errors the previous ones */
        tmp = g->em_n_j;
        g->em_n_j = g->em_n_1_j;
        g->em_n_1_j = tmp;
    }
}
