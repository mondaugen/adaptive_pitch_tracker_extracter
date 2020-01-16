#ifndef PVOC_WINDOWS_H
#define PVOC_WINDOWS_H 

int
get_pvoc_window(float *dest,
           const char *window_type,
           unsigned int length);

unsigned int
pvoc_calc_n_blocks(unsigned int signal_length,
                   unsigned int hop_size,
                   unsigned int window_size);

#endif /* PVOC_WINDOWS_H */
