#ifndef FLOAT_BUF_H
#define FLOAT_BUF_H 
 
struct rngbuf_f32;

void
rngbuf_f32_free(struct rngbuf_f32 *fb);

struct rngbuf_f32 *
rngbuf_f32_new(unsigned int size);

struct rngbuf_f32_where_val { unsigned int n; float f; }; 

int
rngbuf_f32_push_copy(struct rngbuf_f32 *fb, unsigned int n, const float *values);

int
rngbuf_f32_lookup(struct rngbuf_f32 *fb, unsigned int n, float *dest);

/*
process a region starting at start that extends for length.
process may be called multiple times but will always be called in order (say if
the underlying datastructure is segmented).
*/
int
rngbuf_f32_process_region(
    struct rngbuf_f32 *fb,
    unsigned int start,
    unsigned int length,
    void (*process)(
        float *seg,
        unsigned int len,
        void *aux),
    void *aux);

int
rngbuf_f32_shift_in(
    struct rngbuf_f32 *fb,
    const float *values,
    unsigned int nvalues);

int
rngbuf_f32_advance_head(struct rngbuf_f32 *rb,
                        unsigned int n);

/*
Copy values out of the float buffer from start and extending for length into
dest.
*/
int
rngbuf_f32_memcpy(
    struct rngbuf_f32 *fb,
    unsigned int start,
    unsigned int length,
    float *dest);

/*
Make an array of struct rngbuf_f32_where_val.
The value is included in the array if chk returns non-zero when called on the
value.
fun is then called on the array of values.
*/
int
rngbuf_f32_where_values(
    struct rngbuf_f32 *fb,
    unsigned int start,
    unsigned int length,
    int (*chk)(float val, void *aux),
    void (*fun)(struct rngbuf_f32_where_val *v,
                unsigned int nvals,
                void *aux),
    void *aux);

#endif /* FLOAT_BUF_H */
