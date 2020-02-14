#ifndef ALLOC_MMAP_H
#define ALLOC_MMAP_H 

struct alloc_mmap_aux {
    const char *path;
    int fn;
    char *s;
    unsigned int capacity;
};

void alloc_mmap_free(struct alloc_mmap_aux *aux);

void *alloc_mmap(unsigned int capacity,
                 void *_aux);


#endif /* ALLOC_MMAP_H */
