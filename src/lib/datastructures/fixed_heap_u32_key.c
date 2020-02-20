/* Heap of unsigned ints, where the item is just the key */

#include "fixed_heap.h"
#include "fixed_heap_u32_key.h"

static int heap_item_u32_key_min_cmp(void *a_, void *b_, void *aux)
{
    unsigned int *a = a_,
                 *b = b_;
    if (*a > *b) { return 1; }
    return 0;
}

static int heap_item_u32_key_max_cmp(void *a_, void *b_, void *aux)
{
    unsigned int *a = a_,
                 *b = b_;
    if (*b > *a) { return 1; }
    return 0;
}

struct fixed_heap *
fixed_heap_u32_key_new(struct fixed_heap_u32_key_init *init)
{
    struct fixed_heap_init fh_init = {
        .max_n_items = init->max_n_items,
        .item_size = sizeof(unsigned int),
        .cmp = init->max_heap ? heap_item_u32_key_max_cmp :
            heap_item_u32_key_min_cmp,
    };
    return fixed_heap_new(&fh_init);
}
