#include <assert.h>
#include <stdio.h>
#include "test_common.h"
#include "attack_finder.h"

/*
input file should be mono
*/

#define IN_FILE_PATH "/tmp/attacks_from_sd_test_in.f32"
#define OUT_FILE_PATH_BEG "/tmp/attacks_from_sd_test_out_beg.u32"
#define OUT_FILE_PATH_END "/tmp/attacks_from_sd_test_out_end.u32"

int main (void)
{
    float *x = file_to_array(
        IN_FILE_PATH,
        0);
    unsigned int length = get_file_length_path(IN_FILE_PATH)/sizeof(float),
                *beginnings = NULL,
                *ends = NULL;
    struct attacks_from_spec_diff_finder_args args = {
        /* these settings worked well for 16Khz sampled files */
        .H = 256,
        .W = 1024,
        .window_type = "hann",
        .smoothing = 1,
        .lmax_filt_rate = 16000,
        .ng_th = -60,
        .sd_th = 0
    };
    struct attacks_from_spec_diff_finder *finder =
        attacks_from_spec_diff_finder_new(&args);
    assert(finder);
    struct attacks_from_spec_diff_result *ret =
        attacks_from_spec_diff_finder_compute(
        finder,
        x,
        length,
        1);
    assert(ret);
    beginnings = attacks_result_extract_beginnings(ret);
    assert(beginnings);
    ends = attacks_result_extract_ends(ret);
    assert(ends);
    array_to_file(
        OUT_FILE_PATH_BEG,
        beginnings,
        ret->n_attack_time_pairs*sizeof(unsigned int));
    array_to_file(
        OUT_FILE_PATH_END,
        beginnings,
        ret->n_attack_time_pairs*sizeof(unsigned int));
    return 0;
}
