#include <stdio.h>
#include "pvoc_windows.h"

int main (void)
{
    unsigned int H1=100,
                 W1=500,
                 // 0 blocks should fit
                 N=499;

    puts("true 0");
    printf("%d\n",pvoc_calc_n_blocks(N,H1,W1));

    // 1 block should fit
    N=500;
    puts("true 1");
    printf("%d\n",pvoc_calc_n_blocks(N,H1,W1));

    // 1 block should fit
    N=599;
    puts("true 1");
    printf("%d\n",pvoc_calc_n_blocks(N,H1,W1));

    // 2 blocks should fit
    N=600;
    puts("true 2");
    printf("%d\n",pvoc_calc_n_blocks(N,H1,W1));

    // 2 blocks should fit
    N=100;
    puts("true 0");
    printf("%d\n",pvoc_calc_n_blocks(N,H1,W1));

    return 0;
}
