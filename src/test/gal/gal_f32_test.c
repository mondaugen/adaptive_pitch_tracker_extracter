#include <stdio.h>
#include "gal_f32.h"

int main (void)
{

    unsigned int P = 3,
                 N = 10;

    struct gal_f32_init init ={
        .P = P
    };

    float x[] = {9,8,7,6,5,4,3,2,1,0},
          ep[10],
          em[10],
          R[N*P],
          mu[] = {0.99, 0.99, 0.99};

    struct gal_f32_proc proc = {
        .x_in = x,
        .ep = ep,
        .em = em,
        .R  = R,
        .mu = mu,
        .opt = 0,
        .N = 10
    };

    struct gal_f32 *gal = gal_f32_new(&init);
    if (!gal) { fprintf(stderr,"Error initializing gal\n"); }
    gal_f32_process(gal,&proc);
    gal_f32_free(gal);
    return 0;
}
