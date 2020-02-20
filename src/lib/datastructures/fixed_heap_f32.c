#include "fixed_heap_f32.h"

static int heap_item_f32_min_cmp(void *a_, void *b_, void *aux)
{
    struct fixed_heap_item_f32 *a = a_,
                         *b = b_;
    if (a->index > b->index) { return 1; }
    return 0;
}

static int heap_item_f32_max_cmp(void *a_, void *b_, void *aux)
{
    struct fixed_heap_item_f32 *a = a_,
                               *b = b_;
    if (b->index > a->index) { return 1; }
    return 0;
}

struct fixed_heap *
fixed_heap_f32_new(struct fixed_heap_f32_init *init)
{
    struct fixed_heap_init fh_init = {
        .max_n_items = init->max_n_items,
        .item_size = sizeof(struct fixed_heap_item_f32),
        .cmp = init->max_heap ? heap_item_f32_max_cmp : heap_item_f32_min_cmp,
    };
    return fixed_heap_new(&fh_init);
}
