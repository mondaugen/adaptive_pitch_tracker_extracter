#include <stdio.h>
#include <string.h>
#include "ola_f32.h"
#include "test_common.h"

/* Implementation of ola's add */
void ola_f32_add(float *a, const float *b, unsigned int length)
{
    while (length--) {
        *a++ += *b++;
    }
}

int main (void)
{
    float w1[] = {1, 2, 3, 4, 5},
          w2[] =       {5, 6, 7, 8, 9},
          w3[] =             {9, 8, 7, 6, 5},
          w4[] =                   {4, 3, 2, 1,-1},
          w5[] =                        {-2,-3,-4,-5,  -6},
          w6[] =                              {-7,-8,  -9, -8, -7},
          w7[] =                                    {  -6, -5, -4, -3,-2},
          w8[] =                                            {  -1,  9, 8, 7, 6},
          w9[] =                                                    {  5, 4, 3, 2, 1},
          ret[]= {1, 2, 8,10,21,16,20,9, 5,-2,-12,-13,-21,-13,-12,  6,11,11},
          res[sizeof(ret)/sizeof(float)],
          *pres = res,
          *wins[] = {w1,w2,w3,w4,w5,w6,w7,w8,w9},
          **pwins = wins;
    struct ola_f32_init_t config = {
      .sum_in_length = 5,
      .shift_out_length = 2,
    };
    struct ola_f32_t * olaf32 = ola_f32_new(&config);
    int err = 0;
    unsigned int len = sizeof(ret)/sizeof(float), n;
    if (!olaf32) { err = 1; goto fail; }
    for (n = 0; n < len; n += config.shift_out_length) {
        const float *frame = ola_f32_sum_in_and_shift_out(olaf32,*pwins++);
        memcpy(res + n, frame, sizeof(float) * config.shift_out_length);
    }
    if (!CHK_EQ(ret,res,len)) { err = 2; goto fail; }
    printf("correct\n");
    PRINT_ALL(ret,len);
    printf("computed\n");
    PRINT_ALL(res,len);
fail:
    if (olaf32) { ola_f32_free(olaf32); }
    return err;
}



