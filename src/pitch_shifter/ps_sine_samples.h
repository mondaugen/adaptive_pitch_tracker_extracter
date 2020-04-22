#ifndef PS_SINE_SAMPLES_H
#define PS_SINE_SAMPLES_H 
#include "pitch_shifter.h"

/*
NOTE: For testing purposes only because ps_sine_samples_config allocates
memory (which you can't then free)
*/
int
ps_sine_samples_config(float f, struct pitch_shifter_config *config);

#endif /* PS_SINE_SAMPLES_H */
