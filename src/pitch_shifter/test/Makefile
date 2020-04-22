INCLUDE=-I ../../lib/dsp_math/ -I ../../lib/datastructures/ -I ../ -I ../../test/
CPPFLAGS=$(INCLUDE)
CFLAGS=$(INCLUDE) -L ../../lib/dsp_math/ -L ../../lib/datastructures/ -g
deps = ../pitch_shifter.o \
       ../ps_cubic.o \
	   ../ps_sine_samples.o
sources = $(wildcard *.c)

.PHONY=build_libdsp_math

bin/simple_sine_ps: $(sources:.c=.o) build_libdsp_math build_deps build_libdatastructures
	$(CC) $(CFLAGS) $(filter %.o,$^) $(deps) -o $@ -ldsp_math -lm -ldatastructures

include ../../common.mk

include $(sources:.c=.d)

build_deps:
	$(MAKE) -C ../

build_libdsp_math:
	$(MAKE) -C ../../lib/dsp_math/

build_libdatastructures:
	$(MAKE) -C ../../lib/datastructures/
