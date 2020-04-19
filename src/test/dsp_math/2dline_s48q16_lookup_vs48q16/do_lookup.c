#include <stdlib.h>
#include <stdio.h>
#include "test_common.h"
#include "dsp_math.h"

int main (void)
{
    const char input_path[] = "/tmp/dspm_2dline_s48q16_lookup_vs48q16_input.s48q16",
               output_path[] = "/tmp/dspm_2dline_s48q16_lookup_vs48q16_output.s48q16";
    s48q16 x0 = atol(getenv_default("X0","-10"))<<16,
           y0 = atol(getenv_default("Y0","10"))<<16,
           x1 = atol(getenv_default("X1","10"))<<16,
           y1 = atol(getenv_default("Y1","20"))<<16;
    struct dspm_2dline_s48q16 line = dspm_2dline_s48q16_points(x0,y0,x1,y1);
    s48q16 *input = file_to_array(input_path,0);
    long input_length = get_file_length_path(input_path)/sizeof(s48q16);
    dspm_2dline_s48q16_lookup_vs48q16(&line,
                                      input,
                                      (uint32_t)input_length);
    array_to_file(output_path, input, input_length*sizeof(s48q16));
    return 0;
}
