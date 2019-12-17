#include "adsr_envelopes.h"

#undef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))

void
adsr_iir_1st_order_filter(struct adsr_iir_1st_order_filter_args *args)
{
    float yn_1 = *args->yn_1,
            *y = args->y;
    const float *x = args->x,
                 a = args->a;
    unsigned int N = args->N;
    while (N--) {
        *y = *x + a * yn_1;
        yn_1 = *y;
        y++;
        x++;
    }
    *args->yn_1 = yn_1;
}

/* For each value in states equal to state, put 1 in the result, else 0. */
void
adsr_state_eq(struct adsr_state_eq_args *args)
{
    float *result = args->result;
    unsigned int N = args->N;
    const enum adsr_state *states = args->states,
                           state = args->state;
    while (N--) {
        *result = *states == state ? 1. : 0;
        result++;
        states++;
    }
}

void
adsr_subtract_from_const(float *y, float a, const float *b, unsigned int N)
{
    while (N--) {
        *y++ = a - *b++;
    }
}

void
adsr_float_multiply(float *a, const float *b, unsigned int N)
{
    while (N--) {
        *a++ *= *b++;
    }
}

void
adsr_float_add(float *a, const float *b, unsigned int N)
{
    while (N--) {
        *a++ += *b++;
    }
}

void
adsr_extract_trigger(struct adsr_extract_trigger_args *args)
{
    const float *gate = args->gate;
    float last_state = args->last_state,
          *trigger = args->trigger;
    unsigned int N = args->N;
    while (N--) {
        *trigger = *gate - last_state;
        *trigger = MAX(*trigger,0);
        last_state = *gate;
        trigger++;
        gate++;
    }
    args->last_state = last_state;
    adsr_float_multiply(args->trigger,args->sustain_levels,args->N);
}

void
adsr_ramp_smooth(struct adsr_ramp_smooth_args *args)
{
    unsigned int N = args->N, N_mask = args->N - 1;
    float *attack_ramp = args->attack_ramp;
    while (N--) {
        float n = *attack_ramp * args->table_N;
        n = (n < 0) ? 0 : n;
        n = (n > (args->table_N - 1)) ? args->table_N - 1 : n;
        unsigned int floor = (unsigned int)n;
        /* n becomes frac */
        n -= floor;
        /* we do the wrapping because even when frac is 0 (because n was >=
        table_N-1), if one past the end of table is NaN then 0 * NaN is still
        NaN. */
        *attack_ramp = args->table[floor] + n * args->table[(floor+1)&N_mask];
        attack_ramp++;
    }
}
