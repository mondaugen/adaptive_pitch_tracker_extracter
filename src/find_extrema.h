#ifndef FIND_EXTREMA_H
#define FIND_EXTREMA_H 

enum local_max_type {
    local_max_type_right,
    local_max_type_left,
    local_max_type_none,
    local_max_type_both,
};

unsigned int *
local_max_f32(const float *x, 
              unsigned int length, 
              unsigned int *n_max,  
              enum local_max_type type);

unsigned int *
discount_local_max_f32(
    const float *x,
    unsigned int length, 
    unsigned int *n_max, 
    enum local_max_type type, 
    float rate, 
    float min_thresh, 
    float *threshold);

#endif /* FIND_EXTREMA_H */
