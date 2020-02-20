/* A heap with fixed item size */

#include <stdlib.h>
#include <string.h>
#include "fixed_heap.h"

/* assumes item1 and item2 not NULL */
static inline void
item_swap(char *item1, char *item2, unsigned int item_size)
{
    char tmp[item_size];
    memcpy(tmp,item1,item_size);
    memcpy(item1,item2,item_size);
    memcpy(item2,tmp,item_size);
}

static inline unsigned int
parent_index(unsigned int n) { return n & 1 ? ((n+1)>>1) - 1 : (n>>1) - 1; }

static inline unsigned int
left_index(unsigned int n) { return ((n+1)<<1)-1; }

static inline unsigned int
right_index(unsigned int n) { return (n+1)<<1; }

static inline void
heapify(struct fixed_heap *h, unsigned int n)
{
    unsigned int li = left_index(n),
                 ri = right_index(n),
                 ci;
    if (li >= h->cur_n_items) { return; }
    if (ri >= h->cur_n_items) { 
        if (h->cmp(h->items + h->item_size*n, h->items + h->item_size*li,
            h->cmp_aux)) {
            item_swap(h->items + h->item_size*n, h->items + h->item_size*li,
                h->item_size);
        }
        return;
    }
    ci = li;
    if (h->cmp(h->items + h->item_size*li, h->items + h->item_size*ri,
        h->cmp_aux)) {
        ci = ri;
    }
    if (h->cmp(h->items + h->item_size*n, h->items + h->item_size*ci,
        h->cmp_aux)) {
        item_swap(h->items + h->item_size*n, h->items + h->item_size*ci,
            h->item_size);
    }
}

static inline void
float_up(struct fixed_heap *h, unsigned int n)
{
    heapify(h,n);
    if (n != 0) { float_up(h,parent_index(n)); }
}

static inline void
float_down(struct fixed_heap *h, unsigned int n)
{
    if (n >= h->cur_n_items) { return; }
    heapify(h,n);
    float_down(h,left_index(n));
    float_down(h,right_index(n));
}

void
fixed_heap_free(struct fixed_heap *h)
{
    if (h->items) { free(h->items); }
    free(h);
}

struct fixed_heap *
fixed_heap_new(struct fixed_heap_init *init)
{
    struct fixed_heap *ret = calloc(1,sizeof(struct fixed_heap));
    if (!ret) { goto fail; }
    ret->items = calloc(1,init->item_size*init->max_n_items);
    if (!ret->items) { goto fail; }
    ret->max_n_items = init->max_n_items;
    ret->item_size = init->item_size;
    ret->cmp = init->cmp;
    ret->cmp_aux = init->cmp_aux;
    return ret;
fail:
    fixed_heap_free(ret);
    return NULL;
}

/* Returns non-zero if heap is full and item wasn't inserted */
int
fixed_heap_insert(struct fixed_heap *h, void *item)
{
    if (h->max_n_items == h->cur_n_items) { return -1; }
    memcpy(h->items + h->cur_n_items*h->item_size, item, h->item_size);
    h->cur_n_items += 1;
    if (h->cur_n_items > 1) {
        float_up(h,parent_index(h->cur_n_items-1));
    }
    return 0;
}

const void *
fixed_heap_access(struct fixed_heap *h, unsigned int n)
{
    if (n >= h->cur_n_items) { return NULL; }
    return h->items + n*h->item_size;
}

void
fixed_heap_remove_top(struct fixed_heap *h)
{
    if (h->cur_n_items == 0) { return; }
    h->cur_n_items -= 1;
    if (h->cur_n_items > 0) {
        memcpy(h->items,h->items+h->cur_n_items*h->item_size,h->item_size);
        float_down(h,0);
    }
}

    
