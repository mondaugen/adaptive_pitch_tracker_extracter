#include <assert.h>
#include <stdio.h>
#include "test_common.h"

/* included like this so we can access static functions */
#include "attack_finder.c"

/*
input file should be mono
*/

#define IN_FILE_PATH "/tmp/sd_test_in.f32"
#define OUT_FILE_PATH "/tmp/sd_test_out.f32"

int main (void)
{
    float *x = file_to_array(
        IN_FILE_PATH,
        0);
    unsigned int length = get_file_length_path(IN_FILE_PATH)/sizeof(float);
    struct spec_diff_finder_init init = {
        .H = 256,
        .W = 1024,
        .init_window = get_normalized_window,
        .init_window_aux = "hann"
    };
    struct spec_diff_finder *sdf = spec_diff_finder_new(&init);
    assert(sdf);
    struct spec_diff_result *sdr = spec_diff_finder_find(sdf,x,length);
    assert(sdr);
    assert(array_to_file(OUT_FILE_PATH,
    sdr->spec_diff,
    sdr->length * sizeof(sdr->spec_diff[0])) == 0);
    return 0;
}

