#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include "test_common.h"
#include "aug_note_trans.h"

#define BIG_NOTE_STATE_ARRAY_LEN 100000000
#define N 1000

int main(void)
{
    float *note_states = calloc(BIG_NOTE_STATE_ARRAY_LEN,sizeof(float));
    unsigned int n;
    if (!note_states) { return -1; }
    srandom(time(NULL));
    for (n = 0; n < BIG_NOTE_STATE_ARRAY_LEN; n++) {
        note_states[n] = (random()/(float)RAND_MAX) > 0.5;
    }
    for (n = 0; n <= (BIG_NOTE_STATE_ARRAY_LEN-N); n += N) {
        float aug_note_trans1[N+1];
        float aug_note_trans2[N+1];
        compute_aug_note_trans_indexing(note_states,aug_note_trans1,N);
        compute_aug_note_trans_pointer(note_states,aug_note_trans2,N);
        assert(CHK_EQ(aug_note_trans1,aug_note_trans2,N+1));
    }
    return 0;
}
    

    
