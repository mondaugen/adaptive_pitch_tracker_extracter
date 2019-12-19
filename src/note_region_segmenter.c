/*
Take arrays start, end, active of length N and compute regions where a
routine needs to make calls to fill these sections of an array. This kind of
scenario occurs when you have multiple notes per block of samples: at the
beginning of each note, you need to make a new call to some routine to fill that
note.
start contains a 1 where notes start, and 0s otherwise. 
end contains a 1 the first sample after a note, 0s otherwise.
active contains 1s where the note is active. This is equal to the cumulative sum
of start - end.
N is the length of each one of these signals.
Returns regions and number of regions in N_regions.
Regions must contain space for N regions.
*/

#include "note_region_segmenter.h"

#undef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))

static inline void
compute_aug_note_trans(
    const float *note_states, 
    /* must contain N+1 values */
    float *aug_note_trans,
    unsigned int N)
{
    float last_val = 0.;
    unsigned int n;
    for (n = 0; n < N; n++) {
        aug_note_trans[n] = note_states[n] - last_val;
        last_val = note_states[n];
    }
    aug_note_trans[n] = -last_val;
}

static inline void
one_if_gt(float *a, const float *b, const float c, unsigned int N)
{ while (N--) { *a++ = *b++ > c; } }

static inline void
mone_if_lt(float *a, const float *b, const float c, unsigned int N)
{ while (N--) { *a++ = -1 * (*b++ < c); } }

static inline void
pluseq_float(float *a, const float *b, unsigned int N)
{ while (N--) { *a++ += *b++; } }

static inline void
minuseq_float(float *a, const float *b, unsigned int N)
{ while (N--) { *a++ -= *b++; } }

static inline void
where_gt_zero(unsigned int *idcs, const float *a, unsigned int N, unsigned int *N_idcs)
{
    *N_idcs=0;
    unsigned int n;
    for (n = 0; n < N; n++) {
        if (*a++ > 0) {
            *idcs++ = n;
            *N_idcs += 1;
        }
    }
}

static inline void
where_lt_zero(unsigned int *idcs, const float *a, unsigned int N, unsigned int *N_idcs)
{
    *N_idcs=0;
    unsigned int n;
    for (n = 0; n < N; n++) {
        if (*a++ < 0) {
            *idcs++ = n;
            *N_idcs += 1;
        }
    }
}

void
region_segmenter_update(
    const float *start,
    const float *end,
    const float *state,
    unsigned int N,
    struct note_region *regions,
    unsigned int *N_regions)
{
    float aug_note_trans[N+1], aug_note_start[N+1], aug_note_end[N+1];
    unsigned int n, start_idcs[N+1], end_idcs[N+1], N_start_idcs, N_end_idcs;
    compute_aug_note_trans(state,aug_note_trans,N);
    one_if_gt(aug_note_start,aug_note_trans,0,N+1);
    pluseq_float(aug_note_start,start,N);
    one_if_gt(aug_note_start,aug_note_start,1,N+1);
    mone_if_lt(aug_note_end,aug_note_trans,0,N+1);
    minuseq_float(aug_note_end,end,N);
    mone_if_lt(aug_note_end,aug_note_end,-1,N+1);
    where_gt_zero(start_idcs,aug_note_start,N+1,&N_start_idcs);
    where_lt_zero(end_idcs,aug_note_end,N+1,&N_end_idcs);
    *N_regions = 0;
    for (n = 0; n < MIN(N_start_idcs,N_end_idcs); n++) {
        regions->start = start_idcs[n];
        regions->end = end_idcs[n];
        *N_regions += 1;
        regions++;
    }
}
