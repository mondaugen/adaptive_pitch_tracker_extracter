#ifndef FIXED_HEAP_U32_KEY_H
#define FIXED_HEAP_U32_KEY_H 

#include "fixed_heap.h"

struct fixed_heap_u32_key_init {
    unsigned int max_n_items;
    /* if non-zero, is a max heap, otherwise a min heap */
    int max_heap;
};

struct fixed_heap *
fixed_heap_u32_key_new(struct fixed_heap_u32_key_init *init);

#endif /* FIXED_HEAP_U32_KEY_H */
