.PHONY: clean

INCLUDE=-I.
CFLAGS=-g
DEPS=dsp_math.h

src_dir=linux_native

# default FFT library is KISS FFT 
FFT_LIB?=kiss
ifeq ($(FFT_LIB),kiss)
INCLUDE+= -I../kissfft/ -I../kissfft/tools/
ext_objects+=$(src_dir)/kiss_fft.o $(src_dir)/kiss_fftr.o
endif

objects=$(patsubst %.c, %.o, $(wildcard $(src_dir)/*.c))

libdsp_math.a : $(objects) $(ext_objects)
	$(AR) rcs $@ $(objects) $(ext_objects)

# how to build KISS FFT objects
ifeq ($(FFT_LIB),kiss)

$(src_dir)/kiss_fft.o: ../kissfft/kiss_fft.c ../kissfft/kiss_fft.h
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

$(src_dir)/kiss_fftr.o: ../kissfft/tools/kiss_fftr.c ../kissfft/kiss_fft.h
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

endif

$(objects) : $(src_dir)/%.o : $(src_dir)/%.c $(DEPS)
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

clean:
	rm -f libdsp_math.a $(objects)
