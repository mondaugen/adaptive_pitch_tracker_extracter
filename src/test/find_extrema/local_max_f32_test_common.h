#ifndef LOCAL_MAX_F32_TEST_COMMON_H
#define LOCAL_MAX_F32_TEST_COMMON_H 

#include "find_extrema.h"
#include "test_common.h"

static inline int
local_max_test(enum local_max_type type)
{
    const char input_file[] = "/tmp/local_max_f32_input.f32",
               output_file[] = "/tmp/local_max_f32_output.f32";
    unsigned int input_file_len = get_file_length_path(input_file),
                 *lmax = NULL, n_lmax;
    float *x = NULL;
    int ret = 1;
    x = file_to_array(input_file,0);
    if (!x) { goto fail; }
    lmax = local_max_f32(x,input_file_len/sizeof(float), 
                  &n_lmax, type);
    if (!lmax) { goto fail; }
    if (array_to_file(output_file,
        lmax,sizeof(unsigned int) * n_lmax)) {
        goto fail;
    }
    ret = 0;
fail:
    if (x) { free(x); }
    if (lmax) { free(lmax); }
    return ret;
}

#endif /* LOCAL_MAX_F32_TEST_COMMON_H */
