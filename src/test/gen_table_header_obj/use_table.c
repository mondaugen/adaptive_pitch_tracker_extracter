#include <stdio.h>
#include "table.h"

#define N 5

//extern float table;
const float *table_ = table;

int main (void)
{
    int n = 0;
    printf("table_N: %d\n",table_N);
    printf("table address: %p\n",table_);
    printf("contents:\n");
    const float *ptr = table_;
    for (n = 0; n < N; n++) {
        printf("%f ",ptr[n]);
    }
    printf("\n");
    return 0;
}
