#include "dsp_math.h"

static inline void
do_muls(
const float *v0,
const float *v1,
const float *v2,
float c0,
float *x,
uint32_t N)
{
    dspm_mul_vf32_vf32_vf32(x,v0,x,N);
    dspm_mul_vf32_vf32_vf32(x,v1,x,N);
    dspm_mul_vf32_vf32_vf32(x,v2,x,N);
    dspm_mul_vf32_f32(x, c0, N);
}

void
dspm_interp1d4p_vu16q16_vf32_vf32(const u16q16 *xi,
                                  const float *y,
                                  float *yi,
                                  uint32_t N)
{
    uint32_t x0[N];
    float f_2[N],
          f_1[N],
          f[N],
          f1[N],
          tmp[N];
   dspm_floor_vu16q16_vu32(xi,x0,N);
   dspm_sub_vu16q16_vu32_vf32(xi,x0,f,N);
   dspm_add_vf32_f32_vf32(f,-2,f_2,N);
   dspm_add_vf32_f32_vf32(f,-1,f_1,N);
   dspm_add_vf32_f32_vf32(f,1,f1,N);
   /* (f+1)*(f-1)*(f-2)*y[x0]/2 */
   dspm_lookup_vf32_vu32_vf32(y,x0,yi,N);
   do_muls(f_2,f_1,f1,0.5,yi,N);
   /* -(f+1)*f*(f-2)*y[x0+1]/2 */
   dspm_add_vu32_u32(x0,N,1);
   dspm_lookup_vf32_vu32_vf32(y,x0,tmp,N);
   do_muls(f_2,f,f1,-0.5,tmp,N);
   dspm_add_vf32_vf32(yi,tmp,N);
   /* -f*(f-1)*(f-2)*y[x0-1]/6 */
   dspm_add_vu32_u32(x0,N,-2);
   dspm_lookup_vf32_vu32_vf32(y,x0,tmp,N);
   do_muls(f,f_1,f_2,-1./6.,tmp,N);
   dspm_add_vf32_vf32(yi,tmp,N);
   /* (f+1)*f*(f-1)*y[x0+2]/6 */
   dspm_add_vu32_u32(x0,N,3);
   dspm_lookup_vf32_vu32_vf32(y,x0,tmp,N);
   do_muls(f1,f,f_1,1./6.,tmp,N);
   dspm_add_vf32_vf32(yi,tmp,N);
}




