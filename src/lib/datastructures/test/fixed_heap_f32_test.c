#include "test_common.h"
#include "fixed_heap_f32.h"

static inline unsigned int
parent_index(unsigned int n) { return n & 1 ? ((n+1)>>1) - 1 : (n>>1) - 1; }

void print_items_in_heap(struct fixed_heap *heap, unsigned int n_items)
{
    unsigned int n;
    printf("in heap\n");
    for (n = 0; n < n_items; n++) {
        const struct fixed_heap_item_f32 *item = fixed_heap_access(heap,n);
        printf("%u ", item->index);
    }
    printf("\n");
}

void print_heap_valid(struct fixed_heap *heap, unsigned int size)
{
    unsigned int n;
    int ret = 1;
    for (n = 1; n < size; n++) {
        const struct fixed_heap_item_f32 *item = fixed_heap_access(heap,n);
        const struct fixed_heap_item_f32 *parent_item = fixed_heap_access(
            heap,parent_index(n));
        ret &= parent_item->index < item->index;
        if (!ret) {
            printf("failed at %u\n", n);
        }
    }
    if (ret) {
        printf("valid\n");
    } else {
        printf("invalid\n");
    }
}

int main (void)
{
    int ret = 0;
    unsigned int 
                 keys[] = { 7, 13,  4, 13,  1,  7,  2,  4, 10,  4,  5, 11,  9, 12, 12 },
                 //keys[] = { 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 },
                 //keys[] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 },
                 n;
    struct fixed_heap_f32_init init = {
        .max_n_items = 11,
        .max_heap = 0
    };
    struct fixed_heap *heap = fixed_heap_f32_new(&init);
    if (!heap) {
        ret = gen_err_msg("Couldn't allocate heap.");
        goto fail;
    }
    for (n = 0; n < (sizeof(keys)/sizeof(unsigned int)); n++) {
        struct fixed_heap_item_f32 item = {
            .index = keys[n],
            .value = keys[n]
        };
        if (fixed_heap_insert(heap,&item)) {
            break;
        }
        print_items_in_heap(heap,n+1);
    }

        
    printf("popping off heap\n");
    while (1) {
        const struct fixed_heap_item_f32 *item = fixed_heap_access(heap,0);
        //print_items_in_heap(heap,n);
        //print_heap_valid(heap,n);
        n-=1;
        if (!item) { break; }
        printf("%u ", item->index);
        //printf("\n");
        fixed_heap_remove_top(heap);
    }
    printf("\n");
fail:
    if (heap) { fixed_heap_free(heap); }
    print_err_msg(ret);
    return ret;
}
