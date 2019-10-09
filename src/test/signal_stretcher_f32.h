#ifndef SIGNAL_STRETCHER_F32_H
#define SIGNAL_STRETCHER_F32_H 

struct signal_stretcher_f32;

struct signal_stretcher_f32_init {
    float *signal;
    unsigned int signal_length;
    unsigned int window_length;
    const char *window_type;
    unsigned int hop_size;
    float fill;
};

int
signal_stretcher_f32_process(
    const struct signal_stretcher_f32 *ss,
    const int *analysis_points,
    unsigned int n_analysis_points,
    float *output,
    unsigned int output_length);

struct signal_stretcher_f32 *
signal_stretcher_f32_new(struct signal_stretcher_f32_init *ssf32i);

void
signal_stretcher_f32_free(struct signal_stretcher_f32 *ssf32);

#endif /* SIGNAL_STRETCHER_F32_H */
