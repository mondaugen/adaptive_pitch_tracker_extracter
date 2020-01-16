/* Get windows for spectral analysis / synthesis */

#include <math.h>
#include <string.h>

static void
gen_hann_window_f32(float *dest, unsigned int dest_length)
{
    unsigned int n;
    for (n = 0; n < dest_length; n++) {
        dest[n] = 0.5 * (1 - cos(2.*M_PI*n/dest_length));
    }
}

/* dest must have enough memory to hold window of length length */
int
get_pvoc_window(float *dest,
                const char *window_type,
                unsigned int length)
{
    if (strcmp(window_type,"hann") == 0) {
        gen_hann_window_f32(dest,length);
        return 0;
    }
    return -1;
}

/* Calculate the number of blocks that fit within a signal */
unsigned int
pvoc_calc_n_blocks(unsigned int signal_length,
                   unsigned int hop_size,
                   unsigned int window_size)
{
    unsigned int ret = (signal_length < window_size) ? 0 :
        (signal_length - window_size) / hop_size + 1;
    return ret;
}
