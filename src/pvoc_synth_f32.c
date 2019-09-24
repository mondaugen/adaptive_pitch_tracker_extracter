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

struct pvs_f32_init_t {
    /* The signal is localized in time by multiplying this window (assumed 0
    outside of this array of length window_length). */
    float *analysis_window;
    /* The output is multiplied by this synthesis window */
    float *synthesis_window;
    /* The length of the windows */
    unsigned int window_length;
    /* The number of samples between the beginnings of the 2 analysis windows */
    unsigned int hop_size;
    /* A function that is passed the auxilary data structure, and a
    pvs_f32_sample_lookup_t structure, which then fills this with the number of
    samples */
    void (*get_samples)(
        void *aux,
        struct pvs_f32_sample_lookup_t *info);
    /* Auxiliary structure for get_samples */
    void *get_samples_aux;
};

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
}
    
