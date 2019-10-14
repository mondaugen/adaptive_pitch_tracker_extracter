#include "pvoc_synth.h"
#include "pvoc_synth_f32/routines_linux_native.h"
#include <stdio.h>
#include <complex.h>
#include "test_common.h"

/* Test the functions in the func table on f32 and z64 data. */

int main (void)
{
    struct pvs_init_t pvs_init = pvs_f32_init_new();
    float rm[] = {1*2,3*4,5*6,7*8},
          a[]  = {  2,  4,  6,  8},
          b[]  = {1  ,3  ,5  ,7  },
          ra[] = {1+2,3+4,5+6,7+8},
          tmp[4];
    float za[] = {1,2,3,4,5,6},
          zb[] = {8,7,6,5,4,3},
          zc[] = {2,3,4,5,6,7};
    complex float zm[] = {
        ((complex float*)za)[0] * ((complex float*)zb)[0],
        ((complex float*)za)[1] * ((complex float*)zb)[1],
        ((complex float*)za)[2] * ((complex float*)zb)[2],
    };
    complex float zd[] = {
        ((complex float*)zb)[0] / ((complex float*)zc)[0],
        ((complex float*)zb)[1] / ((complex float*)zc)[1],
        ((complex float*)zb)[2] / ((complex float*)zc)[2],
    };
    float rza[] = {  7,  8,  9},
          rzm[] = {2*7, 3*7, 4*8, 5*8, 9*6, 9*7},
          r_abs[] = {
             cabs(((complex float*)rzm)[0]), 
             cabs(((complex float*)rzm)[1]),
             cabs(((complex float*)rzm)[2]),
          };

    pvs_init.func_table->math.real_real_cpymult(
        (struct pvs_real_t *)a,
        (struct pvs_real_t *)b,
        (struct pvs_real_t *)tmp,
        sizeof(a)/sizeof(float));
    printf("real_real_cpymult ? %d\n",CHK_EQ(rm,tmp,sizeof(tmp)/sizeof(float)));
    pvs_init.func_table->math.real_real_mult(
        (struct pvs_real_t *)a,
        (struct pvs_real_t *)b,
        sizeof(a)/sizeof(float));
    printf("real_real_mult ? %d\n",CHK_EQ(rm,a,sizeof(tmp)/sizeof(float)));
    pvs_init.func_table->math.complex_complex_mult(
        (struct pvs_complex_t *)za,
        (struct pvs_complex_t *)zb,
        4);
    printf("complex_complex_mult ? %d\n",CHK_EQ(za,(float*)zm,sizeof(za)/sizeof(float)));
    pvs_init.func_table->math.complex_complex_div(
        (struct pvs_complex_t *)zb,
        (struct pvs_complex_t *)zc,
        4);
    printf("complex_complex_div ? %d\n",CHK_EQ(zb,(float*)zd,sizeof(zb)/sizeof(float)));
    pvs_init.func_table->math.complex_real_mult(
        (struct pvs_complex_t *)zc,
        (struct pvs_real_t *)rza,
        4);
    printf("complex_real_mult ? %d\n",CHK_EQ(rzm,(float*)zc,sizeof(zc)/sizeof(float)));
    pvs_init.func_table->math.complex_abs(
        (struct pvs_complex_t *)rzm,
        (struct pvs_real_t *)tmp,
        4);
    printf("complex_abs ? %d\n",CHK_EQ(r_abs,tmp,sizeof(r_abs)/sizeof(float)));
    return 0;
}
