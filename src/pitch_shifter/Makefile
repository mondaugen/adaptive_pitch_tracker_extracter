INCLUDE=-I../lib/dsp_math/ -I../lib/datastructures/
CPPFLAGS=$(INCLUDE)
CFLAGS=-g
sources = $(wildcard *.c)

.PHONY=all

all: $(sources:.c=.o)

include ../common.mk

include $(sources:.c=.d)
