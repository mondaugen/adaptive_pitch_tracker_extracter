
struct pvs_f32_t {
    struct pvs_init_t config;
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
    if (!init->func_table) { return -6; }
    return 0;
}

