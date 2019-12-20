#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "note_region_segmenter.h"
#include "test_common.h"
#include "adsr_envelopes.h"

#define BLOCK 256
#define SIG_LEN BLOCK*100

void print_region(
    struct note_region *reg,
    FILE *f)
{
    fprintf(f,"%u %u\n",reg->start,reg->end);
}

static inline void
float_array_init_const(float *f, float c, unsigned int N)
{ while (N--) { *f++ = c; } }

static inline void
uint32_array_init_const(unsigned int *u, unsigned int c, unsigned int N)
{ while (N--) { *u++ = c; } }

#define open_write_close(array,path,length)\
    {\
    FILE *fout__ = fopen(path,"w"); \
    fwrite(array,length,1,fout__); \
    fclose(fout__);\
    }

int main(void)
{
    unsigned int n, N_regions, n_reg;
    float *active = file_to_array("/tmp/active.f32",SIG_LEN*sizeof(float)),
          *adsr_envelope = calloc(SIG_LEN,sizeof(float)),
          *start = calloc(SIG_LEN,sizeof(float)),
          *end = calloc(SIG_LEN,sizeof(float)),
          *adsr_active = calloc(SIG_LEN,sizeof(float));
    FILE *fout = fopen("/tmp/regions.txt","w");
    unsigned int attack_duration[BLOCK],
                 decay_duration[BLOCK],
                 release_duration[BLOCK];
    float sustain_level[BLOCK];
    uint32_array_init_const(attack_duration,50,BLOCK);
    uint32_array_init_const(decay_duration,25,BLOCK);
    float_array_init_const(sustain_level,0.5,BLOCK);
    uint32_array_init_const(release_duration,150,BLOCK);
    struct adsr *adsr = adsr_new();
    assert(adsr);
    for (n = 0; n <= (SIG_LEN-BLOCK); n+=BLOCK) {
        adsr_gate_to_adsr_seq_start_end_active_args_alloc(
        adsr_args,
        active+n,
        attack_duration,
        decay_duration,
        sustain_level,
        release_duration,
        adsr_envelope+n,
        BLOCK);
        adsr_gate_to_adsr_seq_start_end_active(adsr,adsr_args);
        adsr_seq_to_env(adsr,adsr_args);
        struct note_region regions[BLOCK];
        region_segmenter_update(
            adsr_args->start,
            adsr_args->end,
            adsr_args->active,
            BLOCK,
            regions,
            &N_regions);
        memcpy(start+n,adsr_args->start,BLOCK*sizeof(float));
        memcpy(end+n,adsr_args->end,BLOCK*sizeof(float));
        memcpy(adsr_active+n,adsr_args->active,BLOCK*sizeof(float));
        for (n_reg = 0; n_reg < N_regions; n_reg++) {
            regions[n_reg].start += n;
            regions[n_reg].end += n;
            print_region(regions+n_reg,fout);
        }
    }
    fclose(fout);
    open_write_close(adsr_envelope,"/tmp/adsr_envelope.f32",SIG_LEN*sizeof(float));
    open_write_close(start,"/tmp/start.f32",SIG_LEN*sizeof(float));
    open_write_close(end,"/tmp/end.f32",SIG_LEN*sizeof(float));
    open_write_close(adsr_active,"/tmp/adsr_active.f32",SIG_LEN*sizeof(float));
    fclose(fout);
    return 0;
}
