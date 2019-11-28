#ifndef FIXED_HEAP_F32_H
#define FIXED_HEAP_F32_H 

#include "fixed_heap.h"

struct fixed_heap_item_f32 {
    unsigned int index;
    float value;
};

struct fixed_heap_f32_init {
    unsigned int max_n_items;
    /* if non-zero, is a max heap, otherwise a min heap */
    int max_heap;
};

struct fixed_heap *
fixed_heap_f32_new(struct fixed_heap_f32_init *init);

#endif /* FIXED_HEAP_F32_H */
