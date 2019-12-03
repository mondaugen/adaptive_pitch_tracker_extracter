#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdint.h>

#define N 10

int fast_floor_log2_f32(float x)
{
    const int emax = -127;
    assert(x>0);
    void *px = &x;
    int exp1 = ((*(uint32_t*)px >> 23) & 0xff) + emax;
    return exp1;
    
}

int main (void)
{
    float x[N];
    int n;
    srand(time(NULL));
    printf("original\n");
    for (n = 0; n < N; n++) {
        x[n] = 10.*rand()/((float)RAND_MAX)+1e-8;
        printf("%f\n",x[n]);
    }
    printf("true floor log2\n");
    for (n = 0; n < N; n++) {
        printf("%f\n",floor(log2(x[n])));
    }
    printf("fast floor log2\n");
    for (n = 0; n < N; n++) {
        printf("%d\n",fast_floor_log2_f32(x[n]));
    }
    return 0;
}


