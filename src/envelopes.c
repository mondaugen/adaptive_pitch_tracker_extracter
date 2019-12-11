#include "envelopes.h"
#include "envelopes_attack_table.h"

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

const struct adsr_init adsr_init_default {
    .decay_min_dB = -60
};

/* all duration values in samples */
struct adsr {
    int attack_duration;
    int decay_duration;
    float sustain_level;
    int release_duration;
    float decay_min_A;
};

