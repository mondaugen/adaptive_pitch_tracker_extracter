#ifndef FIXED_HEAP_H
#define FIXED_HEAP_H 

struct fixed_heap {
    char *items;
    unsigned int cur_n_items;
    unsigned int max_n_items;
    unsigned int item_size;
    int (*cmp)(void *item1, void *item2, void *cmp_aux);
    void *cmp_aux;
};

struct fixed_heap_init {
    unsigned int max_n_items;
    unsigned int item_size;
    /*
    If item1 is the parent of item2 in the heap and cmp(item1,item2,...)
    returns non-zero, then item1 and item2 are swapped.
    */
    int (*cmp)(void *item1, void *item2, void *cmp_aux);
    void *cmp_aux;
};

void fixed_heap_free(struct fixed_heap *h);
struct fixed_heap * fixed_heap_new(struct fixed_heap_init *init);
int fixed_heap_insert(struct fixed_heap *h, void *item);
const void * fixed_heap_access(struct fixed_heap *h, unsigned int n);
void fixed_heap_remove_top(struct fixed_heap *h);

#endif /* FIXED_HEAP_H */
