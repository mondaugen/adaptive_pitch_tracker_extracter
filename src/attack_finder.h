#ifndef ATTACK_FINDER_H
#define ATTACK_FINDER_H 

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
    /* the threshold (> 0) for the spectral difference */
    float sd_th;
};

#define attacks_from_spec_diff_args_default \
(struct attacks_from_spec_diff_finder_args) {\
    .H = 256,\
    .W = 1024,\
    .window_type = "hann",\
    .smoothing = 1,\
    .lmax_filt_rate = 16000,\
    .ng_th=-60,\
    .sd_th=0,\
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
attack_finder_index_mask(
    /* gets freed by this function if it succeeds */
    unsigned int *idcs,
    /* on entry contains length of indcs,
    on exit contains length of masked idcs */
    unsigned int *n_idcs,
    /* all values in gate_changes must be unique 
    and sorted in ascending order */
    const unsigned int *gate_changes,
    /* assumes n_gate >= max(idcs) */
    unsigned int n_gate_changes,
    /* length of signal that gate is based off of */
    unsigned int len_sig);

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

unsigned int *
attacks_result_extract_beginnings(struct attacks_from_spec_diff_result *a);

unsigned int *
attacks_result_extract_ends(struct attacks_from_spec_diff_result *a);

struct attacks_from_spec_diff_finder *
attacks_from_spec_diff_finder_new(
struct attacks_from_spec_diff_finder_args *args);

struct attacks_from_spec_diff_result *
attacks_from_spec_diff_finder_compute(
    struct attacks_from_spec_diff_finder *finder,
    float *x,
    unsigned int length,
    int return_time_pairs);

#endif /* ATTACK_FINDER_H */
