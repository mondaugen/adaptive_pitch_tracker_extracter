#include "dsp_math.h"
#include "test_common.h"
#include <stdio.h>
#include <stdlib.h>

int main (void)
{
    unsigned int N = (unsigned int)get_file_length_path(
    "/tmp/xi.u24q8")/sizeof(float);
    float *y = file_to_array("/tmp/y.f32",0),
          *yi= calloc(N,sizeof(float));
    u24q8 *xi= file_to_array("/tmp/xi.u24q8",0);
    dspm_interp1d4p_vu24q8_vf32_vf32(xi, y, yi, N);
    array_to_file("/tmp/yi_c.f32",yi,N*sizeof(float));
    return 0;
}

