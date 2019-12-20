#include <stdlib.h>
#include <stdio.h>
#include "note_region_segmenter.h"
#include "test_common.h"

#define BLOCK 256
#define SIG_LEN BLOCK*100

void print_region(
    struct note_region *reg,
    FILE *f)
{
    fprintf(f,"%u %u\n",reg->start,reg->end);
}

int main(void)
{
    unsigned int n, N_regions, n_reg;
    float *active = file_to_array("/tmp/active.f32",SIG_LEN*sizeof(float));
    float *start = file_to_array("/tmp/start.f32",SIG_LEN*sizeof(float));
    float *end = file_to_array("/tmp/end.f32",SIG_LEN*sizeof(float));
    FILE *fout = fopen("/tmp/regions.txt","w");
    for (n = 0; n <= (SIG_LEN-BLOCK); n+=BLOCK) {
        struct note_region regions[BLOCK];
        region_segmenter_update(
            start+n,
            end+n,
            active+n,
            BLOCK,
            regions,
            &N_regions);
        for (n_reg = 0; n_reg < N_regions; n_reg++) {
            regions[n_reg].start += n;
            regions[n_reg].end += n;
            print_region(regions+n_reg,fout);
        }
    }
    fclose(fout);
    return 0;
}
