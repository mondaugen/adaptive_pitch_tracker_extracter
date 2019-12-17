#ifndef DSP_MATH_H
#define DSP_MATH_H 

/* Fast floor(log2(x)) */
static inline int
dspm_fast_floor_log2_f32(float x)
{
    unsigned int *ux = (unsigned int*)&x;
    /* Extracts exponent from IEEE 754 floating point number */
    int exp = ((*ux>>23)&0xff) - 127;
    return exp;
}

/* Fast approximate fractional part of log2(x) */
static inline float
dspm_fast_log2_aprox_frac_f32(float x)
{
    unsigned int *ux = (unsigned int*)&x;
    /* (*ux&0x7fffff)**(2**-23) */
    float frac = (*ux&0x7fffff)*1.1920928955078125e-07;
    return frac;
}

    

#endif /* DSP_MATH_H */
