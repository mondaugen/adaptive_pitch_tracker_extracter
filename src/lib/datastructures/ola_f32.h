#ifndef OLA_F32_H
#define OLA_F32_H 

struct ola_f32_t;

struct ola_f32_init_t {
  /* This is usually equal to the window length of the transform */
  unsigned int sum_in_length;
  /* This is usually equal to the hop size or equivalently the audio
     processing block size. This must be a power of 2. */
  unsigned int shift_out_length;
};

/*
Define this so the ola knows how to add. should compute a[n] += b[n] for n in
[0,length)
*/
void ola_f32_add(float *a, const float *b, unsigned int length);
float * ola_f32_sum_in_and_shift_out(struct ola_f32_t *ola, const float *input);
void ola_f32_free(struct ola_f32_t *ola);
struct ola_f32_t * ola_f32_new(struct ola_f32_init_t *config);

#endif /* OLA_F32_H */
