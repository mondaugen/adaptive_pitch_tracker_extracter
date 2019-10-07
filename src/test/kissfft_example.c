/* This shows that KISSFFT scales by N after the forward followed by inverse
transform is applied. */

#include <stdio.h>
#include <complex.h>

#include "kiss_fftr.h"

#define N 8

/* This macro doesn't work */
//    ((t1)a[n]) *= x
#define _ARRAY_MULT(t1,a,x,n) \
    ((t1)a)[n] *= x

int main (void)
{
    kiss_fftr_cfg fftr_cfg_fwd = kiss_fftr_alloc(N,0,NULL,NULL),
                  fftr_cfg_inv = kiss_fftr_alloc(N,1,NULL,NULL);
    kiss_fft_scalar seq[N],
                    seq_res[N];
    float scale = 1./N;
    int n;
    for (n = 0; n < N; n++) {
        seq[n] = n;
    }
    kiss_fft_cpx frq[(N)/2+1];
    kiss_fftr(fftr_cfg_fwd,seq,frq);
    for (n = 0; n < (N/2 + 1); n++) {
        //((complex float*)frq)[n] *= scale;
        _ARRAY_MULT(complex float*,frq,scale,n);
    }
    kiss_fftri(fftr_cfg_inv,frq,seq_res);
    for (n = 0; n < N; n++) {
        printf("%f ", seq_res[n]);
    }
    printf("\n");
    return 0;
}


    
