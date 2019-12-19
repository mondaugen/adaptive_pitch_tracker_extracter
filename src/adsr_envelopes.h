#ifndef ADSR_ENVELOPE_H
#define ADSR_ENVELOPE_H 

enum adsr_state {
    /* outputing 0s */
    adsr_state_Z,
    /* attack state */
    adsr_state_A,
    /* decay state */
    adsr_state_D,
    /* sustain state */
    adsr_state_S,
    /* release state */
    adsr_state_R
};

struct adsr;

struct adsr_iir_1st_order_filter_args {
    /* Output placed here, must have length N. This can point to same memory as x. */
    float *y;
    /* y[n-1] taken from here and put here after function returns */
    float *yn_1;
    /* input taken from here */
    const float *x;
    /* feedback coefficients */
    float *a;
    /* length of x and y */
    unsigned int N;
};

void adsr_iir_1st_order_filter(struct adsr_iir_1st_order_filter_args *args);

struct adsr_gate_to_adsr_seq_start_end_active_args {
    /* All signals have length N */
    /* The signal taking values 0 or 1 that is the gate driving the adsr. */
    const float *gate;
    /* The signal giving the attack time, which is sampled when the attack
    region starts. It will be adjusted to have a minimum value of 1. */
    const unsigned int *attack_duration;
    /* Signal giving the decay time, sampled when decay region starts. Adjusted
    to have minimum value of 1. */
    const unsigned int *decay_duration;
    /* Signal giving the sustain level, sampled when the sustain region starts.
    This can take on any value. */
    const float *sustain_level;
    /* Signal giving the release duration. Adjusted like the other duration signals. */
    const unsigned int *release_duration;
    /*
    This will be filled with the states according to the attack decay sustain
    and release.
    */
    enum adsr_state *adsr_states;
    /* This will have a 1 on the sample where the ADSR envelope starts, and 0 otherwise */
    float *start;
    /* This will have a 1 on the sample where the ADSR envelope ends, and 0 otherwise */
    float *end;
    /* This signal will be 1 while the ADSR envelope is active and 0 otherwise */
    float *active;
    /* This will contain the adsr envelope */
    float *adsr_envelope;
    /* This will contain the attack ramp portion of the ADSR */
    float *attack_ramp;
    /* This will contain the decay portion of the ADSR */
    float *decay_sig;
    /* This will contain the sustain portion of the ADSR */
    float *sustain_sig;
    /* This will contain the release portion of the ADSR */
    float *release_sig;
    /* The number of values in the signals */
    unsigned int N;
};

/* allocate one of these puppies on the stack */
#define adsr_gate_to_adsr_seq_start_end_active_args_alloc(\
    name,\
    gate_,\
    attack_duration_,\
    decay_duration_,\
    sustain_level_,\
    release_duration_,\
    adsr_envelope_,\
    N_)\
    char _ ## name ## _data[\
        sizeof(struct adsr_gate_to_adsr_seq_start_end_active_args)\
        + (adsr_envelope == NULL)*sizeof(float)*N_\
        + 7*sizeof(float)*N_\
        + sizeof(enum adsr_state)*N_];\
    struct adsr_gate_to_adsr_seq_start_end_active_args * name = (void*)_ ## name ## _data;\
    name->gate = gate_;\
    {\
    char *ptr = _ ## name ## _data + \
        sizeof(struct adsr_gate_to_adsr_seq_start_end_active_args);\
    name->attack_duration = attack_duration_;\
    name->decay_duration = decay_duration_;\
    name->sustain_level = sustain_level_;\
    name->release_duration = release_duration_;\
    name->N = N_;\
    if (adsr_envelope_) { name->adsr_envelope = adsr_envelope_; }\
    else { \
        name->adsr_envelope = (void*)ptr; ptr += sizeof(float)*N_;\
    } \
    name->adsr_states = (void*)ptr; ptr += sizeof(enum adsr_state)*N_;\
    name->start = (void*)ptr; ptr += sizeof(float)*N_;\
    name->end = (void*)ptr; ptr += sizeof(float)*N_;\
    name->active = (void*)ptr; ptr += sizeof(float)*N_;\
    name->attack_ramp = (void*)ptr; ptr += sizeof(float)*N_;\
    name->decay_sig = (void*)ptr; ptr += sizeof(float)*N_;\
    name->sustain_sig = (void*)ptr; ptr += sizeof(float)*N_;\
    name->release_sig = (void*)ptr; ptr += sizeof(float)*N_;\
    }


struct adsr_state_eq_args {
    float *result;
    const enum adsr_state *states;
    const enum adsr_state state;
    const unsigned int N;
};

/* For each value in states equal to state, put 1 in the result, else 0. */
void
adsr_state_eq(struct adsr_state_eq_args *args);

struct adsr_extract_trigger_args {
    /* This will contain the release triggers */
    float *trigger;
    /* This takes and then replaces the last gate state */
    float *last_state;
    /* These are the gate states */
    const float *gate;
    /* sustain level sampled on non-zero to zero transition of gate */
    const float *sustain_levels;
    /* length of the signals */
    unsigned int N;
};

void
adsr_extract_trigger(struct adsr_extract_trigger_args *args);

struct adsr_gate_to_ramp_args {
    /* attack gate */
    float *a;
    /* attack durations sampled from here on 0 to non-zero transition of
    attack_duration */
    const unsigned int *attack_duration;
    /* The last gate value of the last block. Will be updated so that it becomes
    the last value for the current block. */
    float *last_gate;
    /* Last gate_to_ramp sum */
    float *last_gtor_cs;
    /* Last attack divisor */
    float *last_attack_div;
    /* length of the arrays */
    unsigned int N;
};

void
adsr_gate_to_ramp(struct adsr_gate_to_ramp_args *args);

struct adsr_ramp_smooth_args {
    /* array of length N containing values between 0 and 1 signifying normalized
    indices at which to look up values in table */
    float *attack_ramp;
    /* length of attack_ramp */
    const unsigned int N;
    /* table where to look up values */
    const float *table;
    /* length of table, this has to be a power of 2, and that condition is not checked! */
    const unsigned int table_N;
};

void
adsr_ramp_smooth(struct adsr_ramp_smooth_args *args);

/* form a += b where a,b vectors */
void 
adsr_float_add(float *a, const float *b, unsigned int N);

/* form a = b + c where a,b,c vectors */
void 
adsr_float_add_out_of_place(float *a, const float *b, const float *c, unsigned int N);

/* form a *= b where a,b vectors */
void
adsr_float_multiply(float *a, const float *b, unsigned int N);

/* form y = a - b where a const, b vector */
void
adsr_subtract_from_const(float *y, float a, const float *b, unsigned int N);

void adsr_free(struct adsr *a);

struct adsr * adsr_new();

void
adsr_gate_to_adsr_seq_start_end_active(
    struct adsr *self,
    struct adsr_gate_to_adsr_seq_start_end_active_args *args);

void adsr_seq_to_env(
    struct adsr *self,
    struct adsr_gate_to_adsr_seq_start_end_active_args *args);

struct adsr_sah_duration_to_coeff_args {
    /* the sample and hold of durations triggers */
    const float *trigger;
    /* The durations */
    const unsigned int *durations;
    /* Once returns, will contain the filter coefficients. Initally this is all
    zeros except for the first value, which contains the last coefficient from
    the last block */
    float *decay_coeffs;
    /* The last decay coefficient */
    float *last_decay_coef;
    /* length of the signals */
    unsigned int N;
    /* A table such that table[0] contains pow(decay_min_A,1) and
    table[decay_coeff_table_length-1] contains
    pow(decay_min_A,1/(pow(2,decay_coeff_table_length-1))) */
    const float *decay_coeff_table;
    /* The length of the table */
    const unsigned int decay_coeff_table_length;
};
void
adsr_sah_duration_to_coeff(struct adsr_sah_duration_to_coeff_args *args);

struct adsr_sah_multiply_sustain_level_args {
    /* where sustain and decay are active */
    float *states;
    /* sustain levels, sampled when states goes from 0 to non-zero */
    const float *sustain_level;
    /* the last sampled sustain level */
    float *last_sustain_level;
    /* the last sustain level state */
    float *last_sustain_level_sah_state;
    /* length of states and sustain_level */
    unsigned int N; 
};
void
adsr_sah_multiply_sustain_level(struct adsr_sah_multiply_sustain_level_args *args);

struct adsr_decay_when_no_gate_args {
    /* Output goes here */
    float *y;
    /* The gate states */
    float *gate;
    /* The past output value */
    float *yn_1;
    /* The decay coefficients */
    float *a;
    /* The length of the signals */
    unsigned int N;
};
void
adsr_decay_when_no_gate(struct adsr_decay_when_no_gate_args *args);

#endif /* ADSR_ENVELOPE_H */
