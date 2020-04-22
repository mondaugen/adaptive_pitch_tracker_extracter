/*
Make a pitch shifter that pitch shifts a synthetic sine tone for testing.
*/

#include "pitch_shifter.h"
#include "ps_cubic.h"
#include "ps_sine_samples.h"
#include "ps_proc_file.h"
#include "test_common.h"

struct pitch_shifter *
simple_sine_ps_new(
    float f,
    float ps_min,
    float ps_max,
    uint32_t B)
{
    struct pitch_shifter_config config = {
        .ps_min = ps_min,
        .ps_max = ps_max,
        .B = B
    };
    if (ps_sine_samples_config(f,&config)) {
        return NULL;
    }
    ps_cubic_config(&config);
    struct pitch_shifter *ret = pitch_shifter_new(&config);
    return ret;
}

int main (void)
{
    /* frequency */
    float F = atof(getenv_default("F","1000.")),
          SR = atof(getenv_default("SR","16000")),
          PS_MIN = atof(getenv_default("PS_MIN","0.5")),
          PS_MAX = atof(getenv_default("PS_MAX","2.0"));
    int B = atoi(getenv_default("B","1024"));
    struct pitch_shifter *ps = simple_sine_ps_new(F/SR,PS_MIN,PS_MAX,(uint32_t)B);
    if (!ps) { return 1; }
    int err = ps_proc_file(ps,"/tmp/ps_rate.u16q16","/tmp/ts_rate.s16q16",
    "/tmp/simple_sine_ps.f32");
    return err;
}
    

