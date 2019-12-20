#include "adsr_envelopes.h"
#include "dsp_math.h"
#include <math.h>
#include <stdio.h>

#undef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#undef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))

void
adsr_iir_1st_order_filter(struct adsr_iir_1st_order_filter_args *args)
{
    float yn_1 = *args->yn_1,
            *y = args->y;
    const float *x = args->x,
                 *a = args->a;
    unsigned int N = args->N;
    while (N--) {
        *y = *x + *a * yn_1;
        yn_1 = *y;
        y++;
        x++;
        a++;
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
    float last_state = *args->last_state,
          *trigger = args->trigger;
    unsigned int N = args->N;
    while (N--) {
        *trigger = *gate - last_state;
        *trigger = MAX(*trigger,0);
        last_state = *gate;
        trigger++;
        gate++;
    }
    *args->last_state = last_state;
    adsr_float_multiply(args->trigger,args->sustain_levels,args->N);
}

void
adsr_ramp_smooth(struct adsr_ramp_smooth_args *args)
{
    unsigned int N = args->N, N_mask = args->table_N - 1;
    float *attack_ramp = args->attack_ramp;
    while (N--) {
        float n = *attack_ramp * args->table_N;
        n = (n < 0) ? 0 : n;
        n = (n > (args->table_N - 1)) ? args->table_N - 1 : n;
        unsigned int flr = (unsigned int)n;
        /* n becomes frac */
        n -= flr;
        /* we do the wrapping because even when frac is 0 (because n was >=
        table_N-1), if one past the end of table is NaN then 0 * NaN is still
        NaN. */
        *attack_ramp = args->table[flr] 
            + n * (args->table[(flr+1)&N_mask] - args->table[flr]);
        attack_ramp++;
    }
}

struct adsr_decay_coeff_lookup_args {
    /* This will be adjusted to be in bounds [1,2**(decay_coeff_table_length-2)] */
    unsigned int decay_time;
    /* A table such that table[0] contains pow(decay_min_A,1) and
    table[decay_coeff_table_length-1] contains
    pow(decay_min_A,1/(pow(2,decay_coeff_table_length-1))) */
    const float *decay_coeff_table;
    /* The length of the table */
    const unsigned int decay_coeff_table_length;
};

static inline float
adsr_decay_coeff_lookup(struct adsr_decay_coeff_lookup_args *args)
{
    args->decay_time = MAX(args->decay_time,1);
    args->decay_time = MIN(args->decay_time,(1 << (args->decay_coeff_table_length - 2)));
    /* table index always non-negative because we forced args->decay_time >= 1 */
    int table_index = dspm_fast_floor_log2_f32(args->decay_time);
    float frac = dspm_fast_log2_aprox_frac_f32(args->decay_time),
          ret = args->decay_coeff_table[table_index] 
                + frac * (args->decay_coeff_table[table_index+1] 
                    - args->decay_coeff_table[table_index]);
    return ret;
}

void
adsr_sah_duration_to_coeff(struct adsr_sah_duration_to_coeff_args *args)
{
    unsigned int n;
    struct adsr_decay_coeff_lookup_args adsr_decay_coeff_lookup_args = {
        .decay_coeff_table = args->decay_coeff_table,
        .decay_coeff_table_length = args->decay_coeff_table_length,
    };
    /* convert to coefficient only where there's a trigger */
    for (n = 0; n < args->N; n++) {
        if (args->trigger[n] != 0) {
            adsr_decay_coeff_lookup_args.decay_time = args->durations[n];
            args->decay_coeffs[n] = adsr_decay_coeff_lookup(&adsr_decay_coeff_lookup_args);
        } else {
            args->decay_coeffs[n] = 0;
        }
    }
    /* now extend these to cover all non-zero values */
    float cur_coef = *args->last_decay_coef;
    for (n = 0; n < args->N; n++) {
        args->decay_coeffs[n] = (args->decay_coeffs[n] == 0) ?
            cur_coef :
            args->decay_coeffs[n]; 
        cur_coef = args->decay_coeffs[n];
    }
    *args->last_decay_coef = cur_coef;
}

static inline void
dur_to_attack_div(
const unsigned int *dur,
const float *trigs,
float *div,
float *last_div,
unsigned int N)
{
    unsigned int n;
    float last_div_ = *last_div;
    for (n = 0; n < N; n++) {
        if (trigs[n] != 0) { last_div_ = 1./MAX(dur[n],1); }
        div[n] = last_div_;
    }
    *last_div = last_div_;
}

void
adsr_gate_to_ramp(struct adsr_gate_to_ramp_args *args)
{
    float last_gate = *args->last_gate,
          last_gtor_cs = *args->last_gtor_cs;
    unsigned int n;
    float dec[args->N], ret[args->N], attack_div[args->N];
    for (n = 0; n < args->N; n++) {
        dec[n] = MIN(last_gate - args->a[n],0);
        last_gate = args->a[n];
    }
    *args->last_gate = last_gate;
    dur_to_attack_div(
        args->attack_duration,
        dec,attack_div,
        args->last_attack_div,
        args->N);
    for (n = 0; n < args->N; n++) {
        ret[n] = last_gtor_cs;
        /*
        Increment last_gtor_cs by 1 when attack gate args->a is non-zero.  Only
        keep past sum if args->a non-zero and a new attack hasn't happened.
        */
        last_gtor_cs = args->a[n] * (1 + last_gtor_cs * (1 + dec[n]));
    }
    *args->last_gtor_cs = last_gtor_cs;
    for (n = 0; n < args->N; n++) {
        args->a[n] = ret[n] * attack_div[n] * args->a[n];
    }
}

void
adsr_sah_multiply_sustain_level(struct adsr_sah_multiply_sustain_level_args *args)
{
    float last_sustain_level = *args->last_sustain_level,
          last_sustain_level_sah_state = *args->last_sustain_level_sah_state,
          diff[args->N];
    unsigned int n;
    for (n = 0; n < args->N; n++) {
        diff[n] = args->states[n] - last_sustain_level_sah_state;
        diff[n] = MAX(diff[n],0);
        last_sustain_level_sah_state = args->states[n];
    }
    *args->last_sustain_level_sah_state = last_sustain_level_sah_state;
    for (n = 0; n < args->N; n++) {
        args->states[n] = diff[n] * args->sustain_level[n]
            + (1 - diff[n]) * last_sustain_level;
        last_sustain_level = args->states[n];
    }
    *args->last_sustain_level = last_sustain_level;
}

void
adsr_decay_when_no_gate(struct adsr_decay_when_no_gate_args *args)
{
    float yn_1 = *args->yn_1;
    unsigned int n;
    for (n = 0; n < args->N; n++) {
        args->y[n] = args->gate[n] == 0 ? args->a[n] * yn_1 : args->gate[n];
        yn_1 = args->y[n];
    }
    *args->yn_1 = yn_1;
}

void 
adsr_float_add_out_of_place(float *a, const float *b, const float *c, unsigned int N)
{
    while (N--) {
        *a++ = *b++ + *c++;
    }
}

void
adsr_extract_start_end_active(struct adsr_extract_start_end_active_args *args)
{
    float last_adsr_start_state = *args->last_adsr_start_state,
          last_adsr_end_state   = *args->last_adsr_end_state,
          trig, state;
    unsigned int n;
    for (n = 0; n < args->N; n++) {
        args->active[n] = args->adsr_states[n] != adsr_state_Z;
    }
    /* extract note starts from beginnings of attack regions */
    for (n = 0; n < args->N; n++) {
        state = args->adsr_states[n] == adsr_state_A;
        trig = state - last_adsr_start_state;
        args->start[n] = trig > 0 ? 1 : 0;
        last_adsr_start_state = state;
    }
    for (n = 0; n < args->N; n++) {
        state = args->adsr_states[n] == adsr_state_R;
        trig = state - last_adsr_end_state;
        args->end[n] = trig < 0 ? 1 : 0;
        last_adsr_end_state = state;
    }
    *args->last_adsr_start_state = last_adsr_start_state;
    *args->last_adsr_end_state = last_adsr_end_state;
}

