#ifndef WINDOWED_LOOKUP_F32_H
#define WINDOWED_LOOKUP_F32_H 

struct windowed_lookup_f32_t;

struct windowed_lookup_f32_init_t {
    /* The size of the lookups */
    unsigned int window_length;
    /* The signal to lookup from. This is not freed when windowed_lookup_f32_t is freed. */
    const float *signal;
    /* The length of this signal */
    unsigned int signal_length;
    /* A function that allows you to fill the out-of-bounds values. This
    function is only called on initialization. If this is NULL, the
    out-of-bounds values are filled with 0s. */
    void (*out_of_bounds_gen)(
        /* Where the values will go */
        float *dest,
        /* How many values to fill */
        unsigned int dest_length,
        /* Auxilary data */
        void *aux);
    /* Auxilary data for out_of_bounds_gen. This is not freed when
    windowed_lookup_f32_t is freed. */    
    void *aux;
};

struct windowed_lookup_f32_access_t {
    /* Section of signal of length window_length */
    const float *signal_section;
    /* Index after access function adjusted it, in the case that the access was
    out of bounds. So the minimum value this could take is -window_length and
    the maxmimum signal_length. */
    const int adjusted_index;
};

void windowed_lookup_f32_free(struct windowed_lookup_f32_t *wl);
struct windowed_lookup_f32_t * windowed_lookup_f32_new(struct windowed_lookup_f32_init_t *config);
void windowed_lookup_f32_out_of_bounds_gen_const(float *dest, unsigned int dest_length, void *aux);
struct windowed_lookup_f32_access_t windowed_lookup_f32_access (struct windowed_lookup_f32_t *wlu, int index);

#endif /* WINDOWED_LOOKUP_F32_H */
