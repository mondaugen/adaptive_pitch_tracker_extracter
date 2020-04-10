#ifndef PITCH_SHIFTER_H
#define PITCH_SHIFTER_H 

struct pitch_shifter_config {
    const float *(*get_samples)(const u24q8 sample_index, void *aux);
    void *get_samples_aux;
    float ps_min;
    float ps_max;
    unsigned int B;
    void (*interpolator)(const u24q8 *xi,
                         const float *y,
                         float *yi,
                         unsigned int N,
                         void *aux);
    void (*interpolator_range)(const u24q8 pos_first,
                               const u24q8 pos_last,
                               int* req_first,
                               int* req_last,
                               void* aux);
    unsigned int (*interpolator_n_points)(unsigned int N, void* aux);
    void *interpolator_aux;
};

#endif /* PITCH_SHIFTER_H */
