/* Casting a function pointer */
#include <stdio.h>

int
my_fun(void *i_, unsigned int len)
{
    int *i = i_, sum = 0;
    while (len--) {
        sum += i[len];
    }
    return sum;
}

int main (void)
{
    int (*fun)(int *, unsigned int) = (int (*)(int *, unsigned int)) my_fun;
    int x[] = {1,2,3,4,5},
        sum = fun(x,5);
    printf("%d\n",sum);
    return 0;
}

/* 1 2 3 4 5 6 7 8 9 10 */
        

