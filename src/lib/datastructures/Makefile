# TODO Add more datastructures as needed

sources = ringbuffer.c rngbuf_f32.c
objects = $(patsubst %.c, %.o , $(sources))

libdatastructures.a : $(objects)
	$(AR) rcs $@ $(objects)

include ../../common.mk

include $(sources:.c=.d)
