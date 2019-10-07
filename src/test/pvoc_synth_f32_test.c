/*
Test src/pvoc_synth.c
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pvoc_synth.h"
#include "pvoc_synth_f32/routines_linux_native.h"
#include "datastructures/windowed_lookup_f32.h"

static void
get_samples (void *aux, struct pvs_real_sample_lookup_t * info)
{
    struct windowed_lookup_f32_t *wl = aux;
    struct windowed_lookup_f32_access_t wla = windowed_lookup_f32_access (
        wl, info->first_sample_index);
    info->samples = (const struct pvs_real_t*)wla.signal_section;
}

static void
gen_scaled_random (float *dest, unsigned int dest_length, void *aux)
{
    float scale = *(float*)aux;
    while (dest_length--) {
        *dest++ = random()/(float)RAND_MAX*scale;
    }
}

/* See if the phase vocoder can sucessfully reproduce the input signal */
static int
chk_identity (void)
{
    const unsigned int signal_length = 8,
                       window_length = 4,
                       hop_size = 2;
    const float analysis_window[] = {0,0.5,1.,0.5},
        /* Includes the scaling for resynthesis */
        //  synthesis_window[] = {0.25,0.25,0.25,0.25}, 
          synthesis_window[] = {0,0.5,1,0.5}, 
        signal[] = {1,2,3,4,5,6,7,8};
        float fill = 1e-6,
            res[signal_length];
    struct windowed_lookup_f32_init_t wli = {
        .window_length = window_length,
        .signal = signal,
        .signal_length = signal_length,
        .out_of_bounds_gen = gen_scaled_random,
        .aux = &fill
    };
    struct windowed_lookup_f32_t *wl = windowed_lookup_f32_new(&wli);
    if (!wl) { goto fail; }
    struct pvs_init_t pvs_init = pvs_f32_init_new();
    pvs_init.user = (struct pvs_user_init_t) {
        .analysis_window = (struct pvs_real_t*) analysis_window,
        .synthesis_window = (struct pvs_real_t*) synthesis_window,
        .window_length = window_length,
        .hop_size = hop_size,
        .get_samples = get_samples,
        .get_samples_aux = wl,
    };
    struct pvs_t * pvs = pvs_new(&pvs_init); 
    if (!pvs) { goto fail; }
    int n;
    for (n = 0; n < signal_length; n += hop_size) {
        const float *cur_frame = (const float*)pvs_process(pvs,n);
        memcpy(res+n,cur_frame,hop_size*sizeof(float));
    }
    for (n = 0; n < signal_length; n++) {
        printf("%f\n",res[n]);
    }
    return 0;
fail:
    return -1;
}

int main (void)
{
    int ret;
    printf("chk_identity\n");
    if ((ret = chk_identity())) {
        printf("failed");
    } else {
        printf("succeeded");
    }
    printf("\n");
    return ret;
}
