#ifndef ENVELOPES_H
#define ENVELOPES_H 

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

struct adsr_init {
    float decay_min_dB;
};

struct adsr;

extern const struct adsr_init_default;

struct adsr_iir_1st_order_filter_args {
    /* Output placed here, must have length N. This can point to same memory as x. */
    float *y;
    /* y[n-1] taken from here and put here after function returns */
    float *yn_1;
    /* input taken from here */
    const float *x;
    /* feedback coefficient */
    float a;
    /* length of x and y */
    unsigned int N;
};

extern void adsr_iir_1st_order_filter(struct adsr_iir_1st_order_filter_args *args);

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

struct adsr_state_eq_args {
    float *result;
    const enum adsr_state *states;
    const enum adsr_state state;
    const unsigned int N;
};

/* For each value in states equal to state, put 1 in the result, else 0. */
extern void
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

extern void
adsr_extract_trigger(struct adsr_extract_trigger_args *args);

#endif /* ENVELOPES_H */
