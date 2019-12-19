void
compute_aug_note_trans_indexing(
    const float *restrict note_states, 
    /* must contain N+1 values */
    float *restrict aug_note_trans,
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

void
compute_aug_note_trans_pointer(
    const float *restrict note_states, 
    /* must contain N+1 values */
    float *restrict aug_note_trans,
    unsigned int N)
{
    float last_val = 0.;
    while (N--) {
        *aug_note_trans++ = *note_states - last_val;
        last_val = *note_states++;
    }
    *aug_note_trans = -last_val;
}


