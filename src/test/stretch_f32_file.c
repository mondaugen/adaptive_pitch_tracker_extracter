/*
Stretch a file containing 32-bit floats.
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "signal_stretcher_f32.h"

#define HUGE_WINDOW_SIZE 10000
/* What to put outside of the signal */
#define FILL 1e-6

const char *getenv_default(
    const char *env_name,
    const char *default_value)
{
    const char *ret = getenv(env_name);
    if (!ret) { return default_value; }
    return ret;
}

/* rewinds to begining of file afterward */
static long get_file_length(
    FILE *f)
{
    fseek(f,0,SEEK_END);
    long ret = ftell(f);
    rewind(f);
    return ret;
}

struct calc_stretched_ret {
    int *analysis_points;
    unsigned int analysis_points_length;
    /* Note that this value is the number of floats not number of bytes */
    long output_file_length;
};

struct calc_stretched_ret *
calculate_stretched_analysis_points(
        long input_file_length,
        unsigned int hop_size,
        float stretch)
{
    long output_file_length = (long)round(input_file_length/stretch);
    unsigned int analysis_points_length = (long)floor(output_file_length/(float)hop_size);
    float analysis_hop_size = hop_size * stretch;
    struct calc_stretched_ret *ret = calloc(1,
        sizeof(struct calc_stretched_ret)+sizeof(int)*analysis_points_length);
    if (!ret) { return NULL; }
    ret->analysis_points = (void*)(ret+1);
    unsigned int n;
    for (n = 0; n < analysis_points_length; n++) {
        ret->analysis_points[n] = (int)round(n*analysis_hop_size);
    }
    ret->analysis_points_length = analysis_points_length;
    ret->output_file_length = output_file_length;
    return ret;
}
    

int main (void)
{
    int ret = 0;
    float *input = NULL, *output = NULL;
    const char *filename = getenv_default("FILENAME","/tmp/in.f32"),
               *output_filename = getenv_default("OUTPUT_FILENAME","/tmp/out.f32"),
               *window_type = getenv_default("WINDOW_TYPE","hann");
    float stretch = atof(getenv_default("S","1."));
    unsigned int window_size = atoi(getenv_default("W","2048")),
                 hop_size = atoi(getenv_default("H","512"));
    FILE *file = NULL, *output_file = NULL;
    struct calc_stretched_ret *calc_stretched_ret = NULL;
    struct signal_stretcher_f32 *signal_stretcher  = NULL;
    if (window_size > HUGE_WINDOW_SIZE) { ret = 5; goto fail; }
    if (hop_size > window_size) { ret = 6; goto fail; }
    if (!filename) { ret = 1; goto fail; }
    file = fopen(filename,"r");
    if (!file) { ret = 2; goto fail; }
    long input_file_length = get_file_length(file);
    input = malloc(input_file_length);
    if (!input) { ret = 3;  goto fail; }
    if ((fread(input,1,input_file_length,file) != input_file_length)) {
        ret = 8; goto fail;
    }
    /* To lengthen a file, we specify analysis points that are closer together
    than the hop size and to shorten, the opposite */
    calc_stretched_ret = calculate_stretched_analysis_points(
        input_file_length/sizeof(float),
        hop_size,
        stretch);
    if (!calc_stretched_ret) { ret = 7; goto fail; }
    output = calloc(calc_stretched_ret->output_file_length,sizeof(float));
    if (!output) { ret = 4; goto fail; }
    struct signal_stretcher_f32_init ssinit = {
        .signal = input,
        .signal_length = input_file_length/sizeof(float),
        .window_length = window_size,
        .window_type = window_type,
        .hop_size = hop_size,
        .fill = FILL
    };
    signal_stretcher = signal_stretcher_f32_new(&ssinit);
    if (!signal_stretcher) { ret = 9; goto fail; }
    if(signal_stretcher_f32_process(
        signal_stretcher,
        calc_stretched_ret->analysis_points,
        calc_stretched_ret->analysis_points_length,
        output,
        calc_stretched_ret->output_file_length)) { ret = 10; goto fail; }
    output_file = fopen(output_filename,"w");
    if (!output_file) { ret = 11; goto fail; }
    if ((fwrite(output,sizeof(float),calc_stretched_ret->output_file_length,output_file)
        != calc_stretched_ret->output_file_length)) {
        ret = 12;
        goto fail;
    }
fail:
    if (file) { fclose(file); }
    if (output_file) { fclose(output_file); }
    if (input) { free(input); }
    if (output) { free(output); }
    if (calc_stretched_ret) { free(calc_stretched_ret); }
    if (signal_stretcher) { signal_stretcher_f32_free(signal_stretcher); }
    return ret;
}
    
    
