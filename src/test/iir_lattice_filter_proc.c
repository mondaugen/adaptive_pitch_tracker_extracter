/* Process f32 files with an IIR lattice filter */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "iir_lattice_filter_f32.h"
#include "test_common.h"

int main (void)
{
    int ret = 0;
    float *in = NULL, *out = NULL,
        *R = NULL;
    FILE *f_in = NULL, *f_out = NULL, *r_file = NULL;
    unsigned int buf_size = atoi(getenv_default("BUF_SIZE","4096")), R_length = 0;
    const char *in_path = getenv_default("IN_PATH","/tmp/in.f32"),
               *out_path = getenv_default("OUT_PATH","/tmp/out.f32"),
               *r_path = getenv_default("R_PATH","/tmp/r.f32");
    const float b0 = atof(getenv_default("B0","1.0"));
    struct iir_lattice_filter_f32 *iirl = NULL;
    f_in = fopen(in_path,"r");
    if (!f_in) { ret = 2; goto fail; }
    f_out = fopen(out_path,"w");
    if (!f_out) { ret = 3; goto fail; }
    r_file = fopen(r_path,"r");
    if (!r_file) { ret = 4; goto fail; }
    R_length = get_file_length(r_file)/sizeof(float);
    R = malloc(sizeof(float)*R_length);
    if (!R) { ret = 5; goto fail; }
    fread(R,sizeof(float),R_length,r_file);
    in = malloc(sizeof(float)*buf_size);
    if (!in) { ret = 6; goto fail; }
    out = malloc(sizeof(float)*buf_size);
    if (!out) { ret = 7; goto fail; }
    struct iir_lattice_filter_f32_init iirli = { .P = R_length };
    struct iir_lattice_filter_f32_proc iirp = {
        .in = in,
        .out = out,
        .N = buf_size,
        .R = R,
        .b0 = b0
    };
    iirl = iir_lattice_filter_f32_new(&iirli);
    if (!iirl) { ret = 1; goto fail; }
    while (!feof(f_in)) {
        memset(in,0,sizeof(float)*buf_size);
        fread(in,sizeof(float),buf_size,f_in);
        iir_lattice_filter_f32_process(iirl,&iirp);
        fwrite(out,sizeof(float),buf_size,f_out);
    }
fail:
    iir_lattice_filter_f32_free(iirl);
    if (in) { free(in); }
    if (out) { free(out); }
    if (R) { free(R); }
    if (f_in) { fclose(f_in); }
    if (f_out) { fclose(f_out); }
    if (r_file) { fclose(r_file); }
    return ret;
}

