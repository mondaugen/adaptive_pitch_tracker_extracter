/*
This object can obtain frames from something providing frames (e.g., an instance
of windowed_lookup_t) and will synthesize recursively as follows:

The last synthesized frame is stored. Then 2 frames H samples apart are obtained
from say a windowed_lookup_t. The output frame is formed by finding the change
in phase from the input frame H samples in the past and the current input frame,
multiplying the output frame by this (to advance its phase) and replacing the
magnitude spectrum of the output frame with that of the current input frame.

In this implementation the analysis and synthesis window must be the same
length.
*/

#include "pvoc_synth_f32.h"
#include "datastructures/rngbuf_f32.h"

struct pvs_f32_t {
    struct pvs_f32_init_t config;
    /* A buffer that output frames are overlapped and added to.
    TODO: It might be faster to use a datastructure specifically made for
    overlapping and adding, becuase the underlying ringbuffer uses a sentinel,
    so accesses will be fragmented. On the otherhand, using the ringbuffer
    allows for non power of 2 window and hop sizes. */
    struct rngbuf_f32 *ola_buffer;
    /* The complex representation of the input at the current time */
    struct pvs_z32_array_t *z_input0;
    /* The complex representation of the input a hopsize ago */
    struct pvs_z32_array_t *z_inputH;
    /* The previous complex representation of the output (updated recursively) */
    struct pvs_z32_array_t *z_outputH;
    /* DFT auxiliary structure */
    struct pvs_f32_dft_t *dft_aux;
    /* Aux processing space, has size config.window_length */
    float *workspace_f32;
};

static int
chk_init_args(struct pvs_f32_dft_init_t *init)
{
    if (!init->analysis_window) { return -1; }
    if (!init->synthesis_window) { return -2; }
    if (init->window_length < 1) { return -3; }
    if (init->hop_size < 1) { return -4; }
    if (!init->get_samples) { return -5; }
    return 0;
}

void
pvs_f32_process(pvs_f32_t *pvs, float *output, int input_time)
{
    /* Get input a hop size ago */
    struct pvs_f32_sample_lookup_t 
        f_inputH = { .first_sample_index = input_time - pvs->config.hop_size,
                     .n_samples = pvs->config.window_length },
        f_input0 = { .first_sample_index = input_time,
                     .n_samples = pvs->config.window_length };
    pvs->get_samples(pvs->get_samples_aux,&f_inputH);
    if (!f_inputH.samples) { return; }
    /* Get current input */
    pvs->get_samples(pvs->get_samples_aux,&f_input0);
    if (!f_input0.samples) { return; }
    /* Perform DFTs */
    pvs_f32_dft_forward(pvs->dft_aux, f_input0.samples, pvs->z_input0);
    pvs_f32_dft_forward(pvs->dft_aux, f_inputH.samples, pvs->z_inputH);
    /* Get magnitude spectrum of 
    /* The previous complex representation of the output (updated recursively) */
    struct pvs_z32_array_t *z_outputH;
}
    
