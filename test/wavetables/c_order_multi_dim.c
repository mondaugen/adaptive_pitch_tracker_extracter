#include <stdio.h>
#define N 10

int main(void)
{
    float x[N*N*N];
    FILE *f = fopen("/tmp/multi_dims.f32","r");
    fread(x,sizeof(float),N*N*N,f);
    int i,j,k;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                printf("%f\n",x[N*N*i+N*j+k]);
            }
        }
    }
    return 0;
}
