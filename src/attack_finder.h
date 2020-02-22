#ifndef ATTACK_FINDER_H
#define ATTACK_FINDER_H 

struct spec_diff_finder_init {
    unsigned int H;
    unsigned int W;
    /* Returns non-zero if it didn't succeed in initializing the window. NOTE:
    From prototype implementation in spectral_difference.py, the window W is
    divided by sum(W) to normalize it, so the init_window should also do this
    normalization to the window it puts in w. */
    int (*init_window)(float *w, void *aux, unsigned int window_length);
    void *init_window_aux;
};

struct attacks_from_spec_diff_finder_args {
    /* hop size */
    unsigned int H;
    /* window size */
    unsigned int W;
    /* window type */
    const char *window_type;
    /*
    smoothing factor s. 0 < s <= 1 1 means no smoothing, 0 would mean totally
    rejecting new values
    */
    float smoothing;
    /*
    max filter discount rate, the time in samples until it reaches 1% of the
    maximum
    */
    unsigned int lmax_filt_rate;
    /* the threshold in dB for the noise gate */
    float ng_th;
};

#define attacks_from_spec_diff_args_default \
(struct attacks_from_spec_diff_finder_args) {\
    .H = 256,\
    .W = 1024,\
    .window_type = "hann",\
    .smoothing = 1,\
    .lmax_filt_rate = 16000,\
    ng_th=-60\
}

/* Used to store time of beginning and end of attack */
struct attack_sample_time_pair {
    unsigned int beginning;
    unsigned int end;
};

struct attacks_from_spec_diff_result {
    struct attack_sample_time_pair *attack_time_pairs;
    unsigned int n_attack_time_pairs;
};

unsigned int *
attack_finder_closest_index_after(
    const unsigned int *filtered,
    /* length of filtered */
    unsigned int n_filtered,
    const unsigned int *find_closest,
    /* length of find_closest */
    unsigned int n_find_closest,
    /* after returning, contains the number of indices in resulting array */
    unsigned int *n_idcs,
    /* if non-zero, filters out indices in find_closest by leaving only those
    coming before indices on the right */ 
    int reverse);

#endif /* ATTACK_FINDER_H */
