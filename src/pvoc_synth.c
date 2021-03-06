#include <stdlib.h>
#include "pvoc_synth.h"

#define PVOC_SMALL_CONST 1e-6
#define MIN(x,y) ((x)<(y)?(x):(y))

struct pvs_t {
    struct pvs_user_init_t config;
    const struct pvs_func_table_t *func_table;
    /* This instant's copies of the analysis and synthesis windows */
    struct pvs_real_t *analysis_window;
    struct pvs_real_t *synthesis_window;
    struct pvs_real_t *init_window;
    /* A buffer that output frames are overlapped and added to. */
    struct pvs_ola_t *ola_buffer;
    /* The complex representation of the input at the current time */
    struct pvs_complex_t *z_input0;
    /* The complex representation of the input a hopsize ago */
    struct pvs_complex_t *z_inputH;
    /* The previous complex representation of the output (updated recursively) */
    struct pvs_complex_t *z_outputH;
    /* DFT auxiliary structure */
    struct pvs_dft_t *dft_aux;
    /* Aux processing space, has size config.window_length */
    struct pvs_real_t *r_workspace;
    /* Space to store output scaling. Has size config.hop_size */
    struct pvs_real_t *output_scaling;
    /* If 0, z_outputH contains garbage because we haven't processed any frames.
    In that case, simply set z_outputH to z_input0 and output z_outputH as if it
    had been computed normally */
    int z_outputH_init;
};

static int
chk_init_args(struct pvs_init_t *init)
{
    struct pvs_user_init_t *uinit = &init->user;
    if (!uinit->analysis_window) { return -1; }
    if (!uinit->synthesis_window) { return -2; }
    if (uinit->window_length < 1) { return -3; }
    if (uinit->hop_size < 1) { return -4; }
    if (!uinit->get_samples) { return -5; }
    /* We assume the functions of func_table are okay, which is reasonable */
    if (!init->func_table) { return -6; }
    return 0;
}

const struct pvs_real_t *
pvs_process(struct pvs_t *pvs, int input_time)
{
    const struct pvs_func_table_t *ftab = pvs->func_table;
    /* Get input a hop size ago */
    struct pvs_real_sample_lookup_t 
        f_inputH = { .first_sample_index = input_time - pvs->config.hop_size,
                     .n_samples = pvs->config.window_length },
        f_input0 = { .first_sample_index = input_time,
                     .n_samples = pvs->config.window_length };
    /* Get current input */
    pvs->config.get_samples(pvs->config.get_samples_aux,&f_input0);
    /* Extract frame and multiply by window */
    ftab->math.real_real_cpymult(
        f_input0.samples,
        pvs->analysis_window,
        pvs->r_workspace,
        pvs->config.window_length);
    if (pvs->z_outputH_init) {
        /* Perform DFT of current input frame */
        ftab->math.dft_forward(
            pvs->dft_aux,
            pvs->r_workspace,
            pvs->z_input0);
        /* Get past input */
        pvs->config.get_samples(pvs->config.get_samples_aux,&f_inputH);
        /* Extract frame and multiply by window */
        ftab->math.real_real_cpymult(
            f_inputH.samples,
            pvs->analysis_window,
            pvs->r_workspace,
            pvs->config.window_length);
        /* Perform DFT of input frame in the past */
        ftab->math.dft_forward(
            pvs->dft_aux,
            pvs->r_workspace,
            pvs->z_inputH);
        /* Find magnitude spectrum of past input frame */
        ftab->math.complex_abs(
            pvs->z_inputH,
            pvs->r_workspace,
            pvs->config.window_length);
        /* Multiply spectrum of current input frame by the magnitude spectrum */
        ftab->math.complex_real_mult(pvs->z_input0,
                                   pvs->r_workspace,
                                   pvs->config.window_length);
        /* Find magnitude spectrum of past output frame */
        ftab->math.complex_abs(
            pvs->z_outputH,
            pvs->r_workspace,
            pvs->config.window_length);
        /* Multiply spectrum of last input frame by the magnitude spectrum */
        ftab->math.complex_real_mult(pvs->z_inputH,
                                   pvs->r_workspace,
                                   pvs->config.window_length);
        /* Find their quotient */
        ftab->math.complex_complex_div(pvs->z_input0,
            pvs->z_inputH,pvs->config.window_length);
        /* Multiply last output spectrum by this quotient */
        ftab->math.complex_complex_mult(
            pvs->z_outputH,pvs->z_input0,pvs->config.window_length);
        /* Inverse fourier transform this output */
        ftab->math.dft_inverse(
            pvs->dft_aux,
            pvs->z_outputH,
            pvs->r_workspace);
        /* Multiply output by synthesis window */
        ftab->math.real_real_mult( pvs->r_workspace,
        pvs->synthesis_window,pvs->config.window_length);
    } else {
        /* Perform DFT of current input frame but put directly in past output frame*/
        ftab->math.dft_forward(
            pvs->dft_aux,
            pvs->r_workspace,
            pvs->z_outputH);
        /*
        Now sum in the first (unwindowed) input frame multiplied by the
        "init_window", which will simulate past inputs having been overlapped and added.
        */
        ftab->math.real_real_cpymult(
            f_input0.samples,
            pvs->init_window,
            pvs->r_workspace,
            pvs->config.window_length);
        pvs->z_outputH_init = 1;
    }
    /* sum into overlap and add buffer shift out samples from overlap and add
    buffer into output */
    struct pvs_real_t *ret = ftab->dstructs.ola_sum_in_and_shift_out(
        pvs->ola_buffer,pvs->r_workspace);
    /* divide-out influence from the windows */
    ftab->math.real_real_mult(ret,pvs->output_scaling,pvs->config.hop_size);
    return ret;
}

void
pvs_free(struct pvs_t *pvs)
{
    if (!pvs) { return; }
    /* pvs should never have been allocated if the config was bad, so we assume
    it is good */
    const struct pvs_func_table_t *ftab = pvs->func_table;
    if (pvs->ola_buffer) { ftab->dstructs.ola_free(pvs->ola_buffer); }
    if (pvs->z_input0) { ftab->dstructs.complex_free(pvs->z_input0); }
    if (pvs->z_inputH) { ftab->dstructs.complex_free(pvs->z_inputH); }
    if (pvs->z_outputH) { ftab->dstructs.complex_free(pvs->z_outputH); }
    if (pvs->dft_aux) { ftab->dstructs.dft_free(pvs->dft_aux); }
    if (pvs->analysis_window) { ftab->dstructs.real_free(pvs->analysis_window); }
    if (pvs->synthesis_window) { ftab->dstructs.real_free(pvs->synthesis_window); }
    if (pvs->init_window) { ftab->dstructs.real_free(pvs->init_window); }
    if (pvs->r_workspace) { ftab->dstructs.real_free(pvs->r_workspace); }
    if (pvs->output_scaling) { ftab->dstructs.real_free(pvs->output_scaling); }
}

struct pvs_t *
pvs_new(struct pvs_init_t *init)
{
    int init_arg_chk_ret;
    if ((init_arg_chk_ret = chk_init_args(init))) {
        return NULL;
    }
    struct pvs_t *ret = calloc(1,sizeof(struct pvs_t));
    if (!ret) { goto fail; }
    ret->config = init->user;
    ret->func_table = init->func_table;
    const struct pvs_func_table_t *ftab = ret->func_table;
    unsigned int W = ret->config.window_length, H = ret->config.hop_size;
    struct pvs_ola_init_t ola_init = {
      .sum_in_length = ret->config.window_length,
      .shift_out_length = ret->config.hop_size
    };
    ret->ola_buffer = ftab->dstructs.ola_alloc(&ola_init);
    if (!ret->ola_buffer) { goto fail; }
    ret->z_input0 = ftab->dstructs.complex_alloc(W);
    if (!ret->z_input0) { goto fail; }
    ret->z_inputH = ftab->dstructs.complex_alloc(W);
    if (!ret->z_inputH) { goto fail; }
    ret->z_outputH = ftab->dstructs.complex_alloc(W);
    if (!ret->z_outputH) { goto fail; }
    struct pvs_dft_init_t dft_init = { .N = W };
    ret->dft_aux = ftab->dstructs.dft_alloc(&dft_init);
    if (!ret->dft_aux) { goto fail; }
    ret->r_workspace = ftab->dstructs.real_alloc(W);
    if (!ret->r_workspace) { goto fail; }
    ret->analysis_window = ftab->dstructs.real_alloc(W);
    if (!ret->analysis_window) { goto fail; }
    ret->synthesis_window = ftab->dstructs.real_alloc(W);
    if (!ret->synthesis_window) { goto fail; }
    ret->init_window = ftab->dstructs.real_alloc(W);
    if (!ret->init_window) { goto fail; }
    ret->output_scaling = ftab->dstructs.real_alloc(H);
    if (!ret->output_scaling) { goto fail; }

    /* Copy analysis and synthesis windows */
    ftab->dstructs.real_memcpy(ret->analysis_window,init->user.analysis_window,W);
    ftab->dstructs.real_memcpy(ret->synthesis_window,init->user.synthesis_window,W);

    /* We don't need the source of the windows any more, so to eliminate
    confusion, we NULL them. */
    ret->config.analysis_window = NULL;
    ret->config.synthesis_window = NULL;

    /* Determine the output scaling before scaling the windows (because we still
    want to cancel out scaling introduced by the DFT */
    unsigned int h;
    for (h = 0; h < W; h += H) {
        struct pvs_real_t *a_win_section = ftab->dstructs.real_offset(
                            ret->analysis_window,h),
                          *s_win_section = ftab->dstructs.real_offset(
                            ret->synthesis_window,h);
        ftab->math.real_real_add_product(
            a_win_section,s_win_section,ret->output_scaling,MIN(H,(W-h)));
    }

    /* Determine the initial window to simulate past overlaps and adds */
    for (h = 0; h < W; h += H) {
        struct pvs_real_t *a_win_section = ftab->dstructs.real_offset(
                            ret->analysis_window,h),
                          *s_win_section = ftab->dstructs.real_offset(
                            ret->synthesis_window,h);
        ftab->math.real_real_add_product(
            a_win_section,s_win_section,ret->init_window,W-h);
    }

    /* Make sure the output_scaling doesn't contain zero or else its reciprocal
    will be bad */
    if (ftab->math.real_contains_zero(ret->output_scaling,H)) { goto fail; }
    ftab->math.real_reciprocal(ret->output_scaling,H);
    /* Scale analysis and synthesis windows */
    struct pvs_dft_window_scale_t dft_window_scale = {
        .dft = ret->dft_aux,
        .analysis_window = ret->analysis_window,
        .synthesis_window = ret->synthesis_window,
        .window_length = W
    };
    ftab->dstructs.dft_window_scale(&dft_window_scale);

    return ret;
fail:
    pvs_free(ret);
    return NULL;
}
