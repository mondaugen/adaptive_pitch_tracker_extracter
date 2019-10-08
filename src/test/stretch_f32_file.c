/*
Stretch a file containing 32-bit floats.
*/
#include <stdlib.h>
#include <math.h>
#include "signal_stretcher_f32.h"

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

int main (void)
{
    int ret = 0;
    float *input = NULL, *output = NULL;
    const char *filename = getenv("FILENAME",NULL),
               *output_filename = getenv("OUTPUT_FILENAME","/tmp/out.f32");
    float stretch = atof(getenv("S","1."));
    FILE *file = NULL, *output_file = NULL;
    if (!filename) { ret = -1; goto fail; }
    file = fopen(filename,"r");
    if (!file) { ret = -2; goto fail; }
    long input_file_length = get_file_length(file),
         output_file_length = (long)round(input_file_length/stretch);;
    input = malloc(sizeof(float)*input_file_length);
    if (!input) { ret = -3;  goto fail; }
    output = malloc(sizeof(float)*output_file_length);
    if (!output) { ret = -4; goto fail; }
    /* TODO add input analysis points etc. */
    struct signal_stretcher_f32_init ssinit = {
        float *signal;
        unsigned int signal_length;
        unsigned int window_length;
        const char *window_type;
        unsigned int hop_size;
        float fill;
    };
fail:
    if (file) { fclose(file); }
    if (output_file) { fclose(output_file); }
    if (input) { free(input); }
    if (output) { free(output); }
    return ret;
}
    
    
