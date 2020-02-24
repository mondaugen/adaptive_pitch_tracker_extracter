#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "fixed_heap_u32_key.h"

#define DEFAULT_MAX 100

int main (void)
{
    struct fixed_heap_u32_key_init hinit = {
        .max_heap = 0
    };
    hinit.max_n_items = 20;
    struct fixed_heap *heap = fixed_heap_u32_key_new(&hinit);
    char *rand_seed_opt = getenv("SEED");
    unsigned int rand_seed = time(NULL), n, val, rmax = DEFAULT_MAX;
    const unsigned int *hval;
    if (rand_seed_opt) {
        assert(sscanf(rand_seed_opt,"%u",&rand_seed)==1);
    }
    srand(rand_seed);
    for (n = 0; n < hinit.max_n_items; n++) {
        val = rand() % DEFAULT_MAX;
        printf("%3u ",val);
        fixed_heap_insert(heap,&val);
    }
    printf("\n");
    while ((hval = fixed_heap_access(heap,0))) {
        val = *hval;
        printf("%3u ",val);
        fixed_heap_remove_top(heap);
    }
    printf("\n");
    fixed_heap_free(heap);
    return 0;
}

