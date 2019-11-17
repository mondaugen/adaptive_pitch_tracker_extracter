/* Process f32 files with a lattice filter */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "lattice_filter_f32.h"
#include "test_common.h"

int main (void)
{
    int ret = 0, lattice_opts = 0;

    float *in = NULL, *out = NULL,
        *R = NULL, *C = NULL;

    FILE *f_in = NULL, *f_out = NULL, *r_file = NULL, *c_file = NULL;

    unsigned int P = atoi(getenv_default("P","0")),
        R_length = 0,
        C_length = 0,
        in_length = 0,
        out_length = 0,
        desired_in_length_for_C = 0,
        desired_in_length_for_R = 0;

    const char *in_path = getenv_default("IN_PATH","/tmp/in.f32"),
               *out_path = getenv_default("OUT_PATH","/tmp/out.f32"),
               *r_path = getenv_default("R_PATH","/tmp/r.f32"),
               *c_path = getenv_default("C_PATH","/tmp/c.f32");

    struct lattice_filter_f32 *latticef = NULL;

    if (P<=0) { ret = gen_err_msg("Bad P"); goto fail; }

    f_in = fopen(in_path,"r");
    if (!f_in) { ret = gen_err_msg("Error opening intput file"); goto fail; }

    f_out = fopen(out_path,"w");
    if (!f_out) { ret = gen_err_msg("Error opening output file"); goto fail; }

    r_file = fopen(r_path,"r");
    if (!r_file) { ret = gen_err_msg("Error opening R file"); goto fail; }

    R_length = get_file_length(r_file)/sizeof(float);
    if (R_length == 0) {
        ret = gen_err_msg("Error getting R file length");
        goto fail;
    } else if (R_length == P) {
        lattice_opts = lattice_filter_f32_NONE;
    } else if ((R_length % P) == 0) {
        lattice_opts = lattice_filter_f32_VARYING_R_C;
    } else {
        ret = gen_err_msg("R length isn't dividable by P");
        goto fail;
    }

    R = malloc(sizeof(float)*R_length);
    if (!R) { ret = gen_err_msg("Error allocating R"); goto fail; }
    fread(R,sizeof(float),R_length,r_file);

    c_file = fopen(c_path,"r");
    if (!c_file) { ret = gen_err_msg("Error opening C file"); goto fail; }
    C_length = get_file_length(c_file)/sizeof(float);
    if (C_length == 0) { ret = gen_err_msg("C length is 0"); goto fail; }
    if ((C_length % (P+1))) {
        ret = gen_err_msg("C_length is not dividable by P+1");
    }
    C = malloc(sizeof(float)*C_length);
    if (!C) { ret = gen_err_msg("Error allocating C"); goto fail; }
    fread(C,sizeof(float),C_length,c_file);

    in_length = get_file_length(f_in)/sizeof(float);
    desired_in_length_for_R = (R_length/P);
    desired_in_length_for_C = (C_length/(P+1));

    if ((lattice_opts & lattice_filter_f32_VARYING_R_C)
        && (in_length != desired_in_length_for_R)) {
        ret = gen_err_msg("Scaled R_length must equal input length");
        goto fail;
    }
    if ((lattice_opts & lattice_filter_f32_VARYING_R_C)
        && in_length != desired_in_length_for_C) {
        ret = gen_err_msg("Scaled C_length must equal input length");
        goto fail;
    }

    in = malloc(sizeof(float)*in_length);
    if (!in) { ret = gen_err_msg("Error allocating input buffer"); goto fail; }

    out = malloc(sizeof(float)*in_length);
    if (!out) { ret = gen_err_msg("Error allocating output buffer"); goto fail; }

    struct lattice_filter_f32_init latticefi = { .P = P };
    struct lattice_filter_f32_proc latticep = {
        .in = in,
        .out = out,
        .N = in_length,
        .R = R,
        .C = C,
        .opts = lattice_opts
    };
    latticef = lattice_filter_f32_new(&latticefi);
    if (!latticef) { ret = gen_err_msg("Error instantiating lattice filter"); goto fail; }
    fread(in,sizeof(float),in_length,f_in);
    lattice_filter_f32_process(latticef,&latticep);
    fwrite(out,sizeof(float),in_length,f_out);
fail:
    lattice_filter_f32_free(latticef);
    if (in) { free(in); }
    if (out) { free(out); }
    if (R) { free(R); }
    if (C) { free(C); }
    if (f_in) { fclose(f_in); }
    if (f_out) { fclose(f_out); }
    if (r_file) { fclose(r_file); }
    if (c_file) { fclose(c_file); }
    print_err_msg(ret);
    return ret;
}

