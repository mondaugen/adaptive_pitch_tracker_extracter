#include <stdlib.h>
#include <string.h>
#include "attack_finder.h"
#include "test_common.h"

#define numelem(a) (sizeof(a)/sizeof(a[0]))

static void *
copy_to_heap(void *p, unsigned int length)
{
    char *ret = malloc(length);
    if (!ret) { return NULL; }
    memcpy(ret,p,length);
    return ret;
}

int main (void)
{
    /*
    signal start high:
    ---___-___
    */
    unsigned int gate_changes_start_high[] = {0,3,6,7},
                 n_gate_changes_start_high = numelem(gate_changes_start_high),
                 idcs_high[] = {0,3,6},
                 n_idcs_high = numelem(idcs_high),
                 correct_idcs_high[] = {0,6},
                 n_correct_idcs_high = numelem(correct_idcs_high),
    /*
    signal start low:
    ___-___---
    */
                 gate_changes_start_low[] = {3,4,7},
                 n_gate_changes_start_low = numelem(gate_changes_start_low),
                 idcs_low[] = {0,3,6,7},
                 n_idcs_low = numelem(idcs_low),
                 correct_idcs_low[] = {3,7},
                 n_correct_idcs_low = numelem(correct_idcs_low),
                 *est_idcs_high = attack_finder_index_mask(
                        copy_to_heap(idcs_high,sizeof(idcs_high)),
                        &n_idcs_high,
                        gate_changes_start_high,
                        n_gate_changes_start_high,
                        10),
                 *est_idcs_low = attack_finder_index_mask(
                        copy_to_heap(idcs_low,sizeof(idcs_low)),
                        &n_idcs_low,
                        gate_changes_start_low,
                        n_gate_changes_start_low,
                        10);
    printf("HIGH PASS? %d\n",
        (n_idcs_high == n_correct_idcs_high) && 
        CHK_EQ(est_idcs_high,correct_idcs_high,n_idcs_high));
    printf("LOW PASS? %d\n",
        (n_idcs_low == n_correct_idcs_low) && 
        CHK_EQ(est_idcs_low,correct_idcs_low,n_idcs_low));
    return 0;
}
