#ifndef AUG_NOTE_TRANS_H
#define AUG_NOTE_TRANS_H 

void
compute_aug_note_trans_indexing(
    const float *restrict note_states, 
    /* must contain N+1 values */
    float *restrict aug_note_trans,
    unsigned int N);

void
compute_aug_note_trans_pointer(
    const float *restrict note_states, 
    /* must contain N+1 values */
    float *restrict aug_note_trans,
    unsigned int N);

#endif /* AUG_NOTE_TRANS_H */
