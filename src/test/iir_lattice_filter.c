/*
Test to see if the lattice filter implementation is correct

The notation used here is consistent with the notation in "Statistical Signal Processing and Modeling" by Monson. See chapter 6.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "test_common.h"

struct iir_lattice_filter_f32 {
    /* Order and number of reflection coefficients */
    unsigned int P;
    /* The P e- values i.e., the backward errors */
    float *e_m;
    /* The P e+ values i.e., the forward errors */
    float *e_p;
    /* The P past e+ values */
    float *e_m_z1;
};

struct iir_lattice_filter_f32_init {
    unsigned int P;
};

void
iir_lattice_filter_f32_free(struct iir_lattice_filter_f32 *l)
{
    if (!l) { return; }
    if (l->e_m) { free(l->e_m); }
    if (l->e_p) { free(l->e_p); }
    if (l->e_m_z1) { free(l->e_m_z1); }
    free(l);
}

struct iir_lattice_filter_f32 *
iir_lattice_filter_f32_new(
    struct iir_lattice_filter_f32_init *init)
{
    if (!init) { return NULL; }
    struct iir_lattice_filter_f32 *ret = calloc(1,sizeof(struct iir_lattice_filter_f32));
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
    iir_lattice_filter_f32_free(ret);
    return NULL;
}

struct iir_lattice_filter_f32_proc {
    const float *in;
    float *out;
    unsigned int N;
    /* P reflection coefficients. The last reflection coefficient is first in
    this array, followed by the second last, etc.*/
    const float *R;
    /* The first and only feedforward coefficient */
    const float b0;
};

void
iir_lattice_filter_f32_process(
    struct iir_lattice_filter_f32 *l,
    struct iir_lattice_filter_f32_proc *p)
{
    float *out = p->out, *e_m, *e_p, *e_m_z1, *tmp;
    const float *in = p->in, *R;
    unsigned int N = p->N, P;
    while (N--) {
        P = l->P;
        e_p = l->e_p + P;
        *e_p = *in * p->b0;
        e_p--;
        R = p->R + P - 1;
        e_m_z1 = l->e_m_z1 + P - 1;
        while (P--) {
            *e_p = *(e_p+1) - *R * *e_m_z1;
            e_p--;
            R--;
            e_m_z1--;
        }
        e_m_z1 = l->e_m_z1;
        e_p = l->e_p;
        *out = *e_p;
        e_m = l->e_m;
        *e_m = *e_p;
        R = p->R;
        P = l->P;
        while (P--) {
            *(e_m + 1) = *e_m_z1 + *R * *e_p;
            e_m++;
            e_m_z1++;
            e_p++;
            R++;
        }
        /* Swap in order to delay */
        tmp = l->e_m;
        l->e_m = l->e_m_z1;
        l->e_m_z1 = tmp;
        in++;
        out++;
    }
}

#define BUF_SIZE 4096

int main (void)
{
    int ret = 0;
    struct iir_lattice_filter_f32_init iirli = { .P = 3 };
    float in[BUF_SIZE], out[BUF_SIZE],
        R[] = { -0.4878, 0.3123, -0.512 };
    struct iir_lattice_filter_f32_proc iirp = {
        .in = in,
        .out = out,
        .N = BUF_SIZE,
        .R = R,
        .b0 = .328
    };
    FILE *f_in = NULL, *f_out = NULL;
    struct iir_lattice_filter_f32 *iirl = iir_lattice_filter_f32_new(&iirli);
    if (!iirl) { ret = 1; goto fail; }
    f_in = fopen("/tmp/in.f32","r");
    if (!f_in) { ret = 2; goto fail; }
    f_out = fopen("/tmp/out.f32","w");
    while (!feof(f_in)) {
        memset(in,0,sizeof(float)*BUF_SIZE);
        fread(in,sizeof(float),BUF_SIZE,f_in);
        iir_lattice_filter_f32_process(iirl,&iirp);
        fwrite(out,sizeof(float),BUF_SIZE,f_out);
    }
fail:
    iir_lattice_filter_f32_free(iirl);
    if (f_in) { fclose(f_in); }
    if (f_out) { fclose(f_out); }
    return ret;
}
