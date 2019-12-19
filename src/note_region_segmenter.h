#ifndef NOTE_REGION_SEGMENTER_H
#define NOTE_REGION_SEGMENTER_H 

struct note_region { unsigned int start; unsigned int end; };

void region_segmenter_update(
    const float *start,
    const float *end,
    const float *state,
    unsigned int N,
    struct note_region *regions,
    unsigned int *N_regions);

#endif /* NOTE_REGION_SEGMENTER_H */
