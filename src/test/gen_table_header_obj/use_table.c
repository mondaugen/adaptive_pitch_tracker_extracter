#include <stdio.h>
#include "table.h"

#define N 5

int main (void)
{
    int n = 0;
    printf("%d\n",table_N);
    for (n = 0; n < N; n++) {
        printf("%d ",table[n]);
    }
    printf("\n");
    return 0;
}
