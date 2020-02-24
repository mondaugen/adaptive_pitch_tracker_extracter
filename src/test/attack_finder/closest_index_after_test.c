#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "attack_finder.h"
#include "test_common.h"

/* Assert x, print error and goto fail if failed */
#define CHK_PRINT_ERR_FAIL(x, msg, ...)\
    if (!x) {\
        fprintf(stderr, msg, ##__VA_ARGS__);\
        fprintf(stderr,"\n");\
        ret = -1;\
        goto fail;\
    }

#define FILTERED_PATH \
    "/tmp/closest_index_after_test_filtered"

#define FIND_CLOSEST_PATH \
    "/tmp/closest_index_after_test_find_closest"

#define CLOSEST_PATH \
    "/tmp/closest_index_after_test_closest"

static int test_find_closest(
    const char *filtered_path,
    const char *find_closest_path,
    const char *closest_path,
    int reverse)
{
    int ret = 0;
    unsigned int *closest = NULL,
                 *filtered = NULL,
                 *find_closest = NULL,
                 n_closest;
    long n_filtered, n_find_closest;
    filtered = file_to_array(filtered_path,0);
    CHK_PRINT_ERR_FAIL(filtered,"loading %s", filtered_path);
    n_filtered = get_file_length_path(
        filtered_path)/sizeof(unsigned int);
    find_closest = file_to_array(find_closest_path,0);
    CHK_PRINT_ERR_FAIL(find_closest,"loading %s", find_closest_path);
    n_find_closest = get_file_length_path(
        find_closest_path)/sizeof(unsigned int);
    closest = attack_finder_closest_index_after(
    filtered,n_filtered,find_closest,n_find_closest,&n_closest,reverse);
    CHK_PRINT_ERR_FAIL(filtered,"calling attack_finder_closest_index_after");
    CHK_PRINT_ERR_FAIL(
        (array_to_file(closest_path,
            closest,n_closest*sizeof(unsigned int)) == 0),
        "saving %s",
        closest_path);
fail:
    if (filtered) { free(filtered); }
    if (find_closest) { free(find_closest); }
    if (closest) { free(closest); }
    return ret;
}

int test_forward()
{
    return test_find_closest(
        FILTERED_PATH "_forward.u32",
        FIND_CLOSEST_PATH "_forward.u32",
        CLOSEST_PATH "_forward.u32",
        0);
}

int test_reverse()
{
    return test_find_closest(
        FILTERED_PATH "_reverse.u32",
        FIND_CLOSEST_PATH "_reverse.u32",
        CLOSEST_PATH "_reverse.u32",
        1);
}

int main (void)
{
    int ret = 0;
    ret = test_forward();
    if (ret) { goto fail; }
    ret = test_reverse();
    if (ret) { goto fail; }
fail:
    return ret;
}
