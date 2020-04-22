/*
Open pitch-shift and time-stretch-rate files and use these to write something to an output file using the pitch-shifter
*/

#include "pitch_shifter.h"
#include "test_common.h"

int
ps_proc_file(
    struct pitch_shifter *ps,
    const char *ps_file_path,
    const char *ts_file_path,
    const char *output_file_path)
{
    u16q16 *ps_rate_sig = NULL;
    s16q16 *ts_rate_sig = NULL;
    float *ret = NULL;
    int err = 0;
    ps_rate_sig = file_to_array(ps_file_path,0);
    if (!ps_rate_sig) { err = -1; goto fail; }
    ts_rate_sig = file_to_array(ts_file_path,0);
    if (!ts_rate_sig) { err = -2; goto fail; }
    long ps_rate_sig_N = get_file_length_path(ps_file_path)/sizeof(u16q16),
         ts_rate_sig_N = get_file_length_path(ts_file_path)/sizeof(s16q16);
    uint32_t rate_sig_N = (uint32_t)(ps_rate_sig_N < ts_rate_sig_N ?
    ps_rate_sig_N : ts_rate_sig_N),
             out_N = 0,
             B = pitch_shifter_B(ps),
             n;
    while (out_N <= (ps_rate_sig_N - B)) { out_N += B; }
    ret = calloc(out_N,sizeof(float));
    if (!ret) { err = -3; goto fail; }
    pitch_shifter_clamp_ps_rate_sig(ps,ps_rate_sig,out_N);
    for (n = 0; n < out_N; n += B) {
        pitch_shifter_process(ps,
                              ps_rate_sig + n,
                              ts_rate_sig + n,
                              ret + n);
    }
    if (array_to_file(output_file_path, ret, out_N*sizeof(float))) { err = -4; goto fail; }
fail:
    if (ps_rate_sig) { free(ps_rate_sig); }
    if (ts_rate_sig) { free(ts_rate_sig); }
    if (ret) { free(ret); }
    return err;
}
