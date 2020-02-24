break attack_finder_closest_index_after
commands
    up
    printf "%p\n", filtered
    eval "print *(unsigned int (*) [%u])filtered", n_filtered
    printf "%p\n", find_closest
    eval "print *(unsigned int (*) [%u])find_closest", n_find_closest
    continue
end

#break fixed_heap_remove_top
#commands
#    printf "%p %u\n", h, *(unsigned int *)h->items
#    continue
#end

break fixed_heap_remove_top
commands
    up
    print n
    printf "%p %u\n", filtered_heap, *(unsigned int *)filtered_heap->items
    printf "%p %u\n", find_closest_heap, *(unsigned int *)find_closest_heap->items
    continue
end

set pagination off
