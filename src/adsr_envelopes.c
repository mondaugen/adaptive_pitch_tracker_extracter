#include "adsr_envelopes.h"
#include "adsr_envelopes_attack_table.h"
#include "adsr_envelopes_decay_coeff_table.h"
#include <stdlib.h>

#undef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))

struct adsr {
    /* all duration values in samples and must be >= 1 */
    unsigned int attack_duration;
    unsigned int decay_duration;
    float sustain_level;
    unsigned int release_duration;
    /* attack count */
    unsigned int attack_n;
    /* decay count */
    unsigned int decay_n;
    /* decay coefficient, calculated from decay_min_A and decay_duration */
    float decay_coeff;
    /* last decay feedback value */
    float decay_yn_1;
    /* release count */
    unsigned int release_n;
    /* release coefficient, calculated from decay_min_A and release_duration */
    float release_coeff;
    /* last release feedback value */
    float release_yn_1;
    /* last adsr state */
    enum adsr_state last_adsr_state;
    /* last gate value */
    float last_gate;
    /* last gate_to_ramp sum */
    float gtor_cs;
    /* last decay state value (0 or 1) */
    float last_decay_state;
    /* last release state value (0 or 1) */
    float last_release_state;
    /* last attack divisor */
    float last_attack_div;
    /* last sampled sustain level */
    float last_sustain_level;
    /* last state for sustain level sampling */
    float last_sustain_level_sah_state;
    /* Last state for differencing filter extracting the start and end of when
    the ADSR is on */
    float last_adsr_gate_state;
};

const struct adsr adsr_default = {
    .attack_duration = 1,
    .decay_duration = 1,
    .sustain_level = 1,
    .release_duration = 1,
    .attack_n = 0,
    .decay_n = 0,
    .decay_yn_1 = 0,
    .release_n = 0,
    .release_yn_1 = 0,
    .last_adsr_state = adsr_state_Z,
    .last_gate = 0,
    .gtor_cs = 0,
    .last_decay_state = 0,
    .last_attack_div = 1,
    .last_sustain_level = 0,
    .last_sustain_level_sah_state = 0,
    .last_adsr_gate_state = 0,
};

void
adsr_free(struct adsr *a)
{
    free(a);
}

struct adsr *
adsr_new(void)
{
    struct adsr *ret = calloc(1,sizeof(struct adsr));
    if (!ret) { goto fail; }
    *ret = adsr_default;
    return ret;
fail:
    if (ret) { adsr_free(ret); }
    return NULL;
}

void
adsr_gate_to_adsr_seq_start_end_active(
    struct adsr *self,
    struct adsr_gate_to_adsr_seq_start_end_active_args *args)
{
    unsigned int n;
    for (n = 0; n < args->N; n++) {
        float g = args->gate[n];
        if (g == 1) {
            if (self->last_adsr_state == adsr_state_Z) {
                self->last_adsr_state = adsr_state_A;
                self->attack_duration = MAX(1,args->attack_duration[n]);
                self->attack_n = 0;
                args->start[n] = 1;
            }
        }
        if (g == 0) {
            if ((self->last_adsr_state != adsr_state_Z) 
                && (self->last_adsr_state != adsr_state_R)) {
                self->last_adsr_state = adsr_state_R;
                self->release_duration = MAX(1,args->release_duration[n]);
                self->release_n = 0;
            }
        }
        args->adsr_states[n]=self->last_adsr_state;
        if (self->last_adsr_state != adsr_state_Z) {
            args->active[n]=1;
        }
        if (self->last_adsr_state == adsr_state_A) {
            self->attack_n += 1;
            if (self->attack_n >= self->attack_duration) {
                self->last_adsr_state = adsr_state_D;
                self->decay_duration = MAX(1,args->decay_duration[n]);
                self->decay_n = 0;
            }
        }
        if (self->last_adsr_state == adsr_state_D) {
            self->decay_n += 1;
            if (self->decay_n >= self->decay_duration) {
                /*
                NOTE the args->sustain_level gets sampled in the
                adsr_seq_to_env function.
                */
                self->last_adsr_state = adsr_state_S;
            }
        }
        if (self->last_adsr_state == adsr_state_R) {
            self->release_n += 1;
            if (self->release_n >= self->release_duration) {
                self->last_adsr_state = adsr_state_Z;
                args->end[n]=1;
            }
        }
    }
}

void adsr_seq_to_env(
    struct adsr *self,
    struct adsr_gate_to_adsr_seq_start_end_active_args *args)
{

    float a[args->N], d[args->N], s[args->N], r[args->N],
         *adsr_gates[] = {a,d,s,r,NULL}, **adsr_gate_ptr = adsr_gates,
         a_ramp[args->N], d_trig[args->N], r_trig[args->N], one_minus_sus_level[args->N],
         decay_coeffs[args->N], sus_level[args->N], ads_states[args->N];
    enum adsr_state adsr_states[] = {adsr_state_A,adsr_state_D,adsr_state_S,adsr_state_R},
         *adsr_state_ptr = adsr_states;
    unsigned int n;
    
    /*
    extract the state sections to signals that are 1 when in that section and
    0 otherwise
    */
    while (*adsr_gate_ptr) {
        struct adsr_state_eq_args adsr_state_eq_args = {
            .result = *adsr_gate_ptr,
            .states = args->adsr_states,
            .state = *adsr_state_ptr,
            .N = args->N,
        };
        adsr_gate_ptr++;
        adsr_state_ptr++;
        adsr_state_eq(&adsr_state_eq_args);
    }

    adsr_float_add_out_of_place(ads_states,a,d,args->N);
    adsr_float_add(ads_states,s,args->N);

    /* Form sustain signal */
    adsr_float_add_out_of_place(sus_level,s,d,args->N);

    /* sample and hold sustain levels */
    struct adsr_sah_multiply_sustain_level_args adsr_sah_multiply_sustain_level_args = {
        .states = sus_level, 
        .sustain_level = args->sustain_level, 
        .last_sustain_level = &self->last_sustain_level, 
        .last_sustain_level_sah_state = &self->last_sustain_level_sah_state,
        .N = args->N, 
    };
    adsr_sah_multiply_sustain_level(&adsr_sah_multiply_sustain_level_args);
    
    /* convert attack gate to ramp */
    struct adsr_gate_to_ramp_args adsr_gate_to_ramp_args = {
        /* attack gate */
        .a = a,
        /* attack durations sampled from here on 0 to non-zero transition of
        attack_duration */
        .attack_duration = args->attack_duration,
        /* The last gate value of the last block */
        .last_gate = &self->last_gate,
        /* The last gate_to_ramp sum */
        .last_gtor_cs = &self->gtor_cs,
        /* Last attack divisor */
        .last_attack_div = &self->last_attack_div,
        /* length of the arrays */
        .N = args->N
    };
    adsr_gate_to_ramp(&adsr_gate_to_ramp_args);

    /* convert attack ramp to one specified by look up table */
    struct adsr_ramp_smooth_args adsr_ramp_smooth_args = {
        .attack_ramp = a,
        .N = args->N,
        /* declared in envelopes_attack_table.h and linked from
        envelopes_attack_table.o */
        .table = adsr_envelopes_attack_table,
        /* declared in envelopes_attack_table.h */
        .table_N = adsr_envelopes_attack_table_length
    };
    adsr_ramp_smooth(&adsr_ramp_smooth_args);

    /* form 1 - sustain level */
    adsr_subtract_from_const(one_minus_sus_level,1,sus_level,args->N);

    /* Convert decay section to impulse marking 0 to non-zero transitions */
    struct adsr_extract_trigger_args adsr_extract_decay_trigger_args = {
        /* This will contain the decay triggers */
        .trigger = d_trig,
        /* This takes and then replaces the last gate state */
        .last_state = &self->last_decay_state,
        .gate = d,
        /* 1 - sustain level sampled on non-zero to zero transition of decay_gate */
        .sustain_levels = one_minus_sus_level,
        .N = args->N,
    };
    adsr_extract_trigger(&adsr_extract_decay_trigger_args);

    /* Sample decay_duration array when d_trig non-zero, convert these to decay
    coefficients */ 
    struct adsr_sah_duration_to_coeff_args decay_dur_to_coeff_args = {
        .trigger = d_trig,
        .durations = args->decay_duration,
        .last_decay_coef = &self->decay_coeff,
        .decay_coeffs = decay_coeffs,
        .N = args->N,
        .decay_coeff_table = adsr_envelopes_decay_coeff_table,
        .decay_coeff_table_length = adsr_envelopes_decay_coeff_table_length,
    };
    adsr_sah_duration_to_coeff(&decay_dur_to_coeff_args);

    /* Filter decay triggers to get decay sections */
    struct adsr_iir_1st_order_filter_args decay_filter_args = {
        .y = d_trig,
        .yn_1 = &self->decay_yn_1,
        .x = d_trig,
        .a = decay_coeffs,
        .N = args->N
    };
    adsr_iir_1st_order_filter(&decay_filter_args);

    /* Multiply the decay filter output by the decay gate */
    adsr_float_multiply(d_trig,d,args->N);

    adsr_float_add(s,d,args->N);
    adsr_float_multiply(s,sus_level,args->N);

    /* Form sum of attack decay and sustain signal */
    adsr_float_add_out_of_place(args->adsr_envelope,a,d_trig,args->N);
    adsr_float_add(args->adsr_envelope,s,args->N);

    /*
    We have to look up release_duration whenever a,d, or s non-zero because the
    next sample could always be the start of the release section, where we'll
    need to know the release duration.
    */ 
    struct adsr_sah_duration_to_coeff_args release_dur_to_coeff_args = {
        .trigger = ads_states,
        .durations = args->release_duration,
        .last_decay_coef = &self->release_coeff,
        .decay_coeffs = decay_coeffs,
        .N = args->N,
        .decay_coeff_table = adsr_envelopes_decay_coeff_table,
        .decay_coeff_table_length = adsr_envelopes_decay_coeff_table_length,
    };
    adsr_sah_duration_to_coeff(&release_dur_to_coeff_args);

    /* Convert release section to impulse marking 0 to non-zero transitions */
    struct adsr_decay_when_no_gate_args adsr_decay_when_no_gate_args = {
        .y = r_trig,
        .gate = args->adsr_envelope,
        .yn_1 = &self->release_yn_1,
        .a = decay_coeffs,
        .N = args->N
    };
    adsr_decay_when_no_gate(&adsr_decay_when_no_gate_args);
    
    /* Multiply the release filter output by the release gate */
    adsr_float_multiply(r_trig,r,args->N);

    /* Add the release section into ads */
    adsr_float_add(args->adsr_envelope,r_trig,args->N);
}
