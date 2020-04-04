#include "dsp_math.h"
#include "test_common.h"
#include <stdio.h>
#include <stdlib.h>

int main (void)
{
    unsigned int N = (unsigned int)get_file_length_path(
    "/tmp/xi.f32")/sizeof(float);
    float *y = file_to_array("/tmp/y.f32",0),
          *xi= file_to_array("/tmp/xi.f32",0),
          *yi= calloc(N,sizeof(float));
    dspm_interp1d4p_vf32_vf32_vf32(xi, y, yi, N);
    array_to_file("/tmp/yi_c.f32",yi,N*sizeof(float));
    return 0;
}

