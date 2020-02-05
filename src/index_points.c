#include <stddef.h>
#include <stdlib.h>
#include "index_points.h"

struct index_points *
index_points_new(unsigned int length)
{
    char *ret = malloc(
    sizeof(struct index_points)+sizeof(unsigned int [length]));
    if (!ret) { return NULL; }
    *(unsigned int *)(ret + offsetof(struct index_points,length)) = length;
    return (struct index_points*)ret;
}

void
index_points_free(struct index_points *x)
{
    free(x);
}
