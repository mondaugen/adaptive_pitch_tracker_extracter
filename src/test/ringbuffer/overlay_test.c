#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include "ringbuffer.h"
#include "alloc_mmap.h"

/*
Code for sending some stuff using a ringbuffer (lock-free queue).
*/

#define RB_CAPACITY 256

const char *strs[] = {
    "hey",
    "we",
    "are",
    "a",
    "happy",
    NULL
};

const char **cur_str = strs;

int running = 1;

static void send_next_string(struct rngbuf *rb)
{
    /* Get length of formatted string */
    int slen = snprintf(NULL,0,"%s\n",*cur_str);
    /* Allocate enough memory */
    char buf[slen+1];
    /* Format it */
    sprintf(buf,"%s\n",*cur_str);
    if ((rngbuf_push_copy(rb,buf,slen) == 0)) {
        cur_str = *(cur_str + 1) == NULL ? strs : cur_str + 1;
    }
}

static void
end_prog(int sn)
{
    running = 0;
}

int main (void)
{
    int ret = 0;
    struct alloc_mmap_aux alloc_mmap_aux = {
        .path = "overlay_test_output_path"
        /* the rest gets filled in by alloc_mmap */
    };
    struct rngbuf *rb = rngbuf_overlay(RB_CAPACITY,
    alloc_mmap,&alloc_mmap_aux);
    if (!rb) { ret = -1; goto fail; }
    signal(SIGINT,end_prog);
    while (running) {
        send_next_string(rb);
        /* wait 10 ms */
        usleep(10000);
    }
fail:
    if (rb) { alloc_mmap_free(&alloc_mmap_aux); }
    return ret;
}

