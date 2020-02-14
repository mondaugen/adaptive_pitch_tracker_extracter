/*
Some slightly hacky routines so you can overlay a ringbuffer onto memory
obtained via mmap. That way you can use the ringbuffer to do inter-process
communication, which is handy because the ringbuffer is lock-free for single
producer and signal consumer threads.
*/

#include <stdlib.h>
#include <stdio.h>
#include <sys/mman.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "alloc_mmap.h"

#define report_and_fail(s) perror(s); goto fail

/*
Call this on the aux you used when you used alloc_mmap to free the memory it
allocated
*/
void alloc_mmap_free(struct alloc_mmap_aux *aux)
{
    if (!aux) { return; }
    if (aux->fn >= 0) { close(aux->fn); }
    if (aux->s) { munmap(aux->s,aux->capacity); }
}

void *alloc_mmap(
    unsigned int capacity,
    void *_aux)
{
    struct alloc_mmap_aux *aux = _aux;
    aux->capacity = capacity;
    aux->fn = -1;
    aux->s = NULL;
    aux->fn = shm_open(aux->path, O_RDWR | O_CREAT, 0644);
    if (aux->fn < 0) { report_and_fail("opening shared memory"); }
    ftruncate(aux->fn,capacity);
    aux->s = mmap(NULL, capacity,
    PROT_WRITE | PROT_READ,
    MAP_SHARED, aux->fn, 0);
    if ((aux->s == MAP_FAILED)) { aux->s = NULL; report_and_fail("mapping memory"); }
    return aux->s;
fail:
    alloc_mmap_free(aux);
    return 0;
}
    
