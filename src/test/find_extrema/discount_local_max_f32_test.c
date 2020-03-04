#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "test_common.h"
#include "find_extrema.h"

#define IN_FILE "/tmp/noi.f32"
#define OUT_FILE "/tmp/noi"
/* Number of samples to reach .01 of maximum */
#define LMFR_N 500
/* threshold above which maximum is accepted */
#define MIN_THRESH 0

int main (void)
{
    const unsigned int len_x = get_file_length_path(IN_FILE)/sizeof(float);
    float *x = file_to_array("/tmp/noi.f32",0),
          *thresh = malloc(len_x*sizeof(float)),
          rate = pow(.01,1./LMFR_N);
    assert(x);
    assert(thresh);
    unsigned int n_max,
                //*lmaxs = local_max_f32(x,len_x,&n_max,local_max_type_right);
                *lmaxs = discount_local_max_f32(
                 x,len_x,&n_max,local_max_type_right,rate,MIN_THRESH,thresh);
    assert(lmaxs);
    array_to_file(OUT_FILE "thresh.f32",thresh,len_x*sizeof(float));
    array_to_file(OUT_FILE "maxs.u32",lmaxs,n_max*sizeof(unsigned int));
    return 0;
}

