src/test/bin/ringbuffer_overlay_test : \
    src/test/ringbuffer/overlay_test.c \
    src/lib/datastructures/ringbuffer.c \
    src/lib/datastructures/ringbuffer.h \
    src/test/alloc_mmap.c \
    src/test/alloc_mmap.h
	cc -g -Isrc -Isrc/test -Isrc/lib/datastructures \
    $(filter %.c,$^) -o $@ -lrt
