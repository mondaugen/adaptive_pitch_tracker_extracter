#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "test_common.h"
#include "adsr_envelopes.h"

#define BLOCK_SIZE 256
int main (void)
{
    long longest_length = 0;
    char *file_paths[]={
        "/tmp/gate.f32",
        "/tmp/attack_duration.u32",
        "/tmp/decay_duration.u32",
        "/tmp/sustain_level.f32",
        "/tmp/release_duration.u32",
        NULL},
        **ptr_file_path=file_paths;
    unsigned int type_sizes[] = {
        sizeof(float),
        sizeof(unsigned int),
        sizeof(unsigned int),
        sizeof(float),
        sizeof(unsigned int)},
       *ptr_type_sizes=type_sizes;
    while (*ptr_file_path) {
        long length = get_file_length_path(*ptr_file_path) / *ptr_type_sizes;
        if (length > longest_length) { longest_length = length; };
        ptr_file_path++;
        ptr_type_sizes++;
    }
    const float        *gate = file_to_array(
                                "/tmp/gate.f32",
                                longest_length*sizeof(float));
    const unsigned int *attack_duration = file_to_array(
                                            "/tmp/attack_duration.u32",
                                            longest_length*sizeof(unsigned int));
    const unsigned int *decay_duration = file_to_array(
                                            "/tmp/decay_duration.u32",
                                            longest_length*sizeof(unsigned int));
    const float        *sustain_level = file_to_array(
                                            "/tmp/sustain_level.f32",
                                            longest_length*sizeof(float));
    const unsigned int *release_duration = file_to_array(
                                            "/tmp/release_duration.u32",
                                            longest_length*sizeof(unsigned int));
    float *adsr_envelope = calloc(longest_length,sizeof(float));
    unsigned long n;
    struct adsr *adsr = adsr_new();
    assert(adsr);
    for (n = 0; n <= (longest_length - BLOCK_SIZE); n += BLOCK_SIZE) {
        adsr_gate_to_adsr_seq_start_end_active_args_alloc(
        adsr_args,
        gate+n,
        attack_duration+n,
        decay_duration+n,
        sustain_level+n,
        release_duration+n,
        adsr_envelope+n,
        BLOCK_SIZE);
        adsr_gate_to_adsr_seq_start_end_active(adsr,adsr_args);
        adsr_seq_to_env(adsr,adsr_args);
    }
    FILE *f = fopen("/tmp/adsr_envelope.f32","w");
    fwrite(adsr_envelope,sizeof(float),longest_length,f);
    fclose(f);
    return 0;
}
