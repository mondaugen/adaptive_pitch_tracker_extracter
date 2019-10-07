src/test/bin/kiss_fft_example : src/test/kissfft_example.c \
                                src/lib/kissfft/kiss_fft.c \
                                src/lib/kissfft/tools/kiss_fftr.c
	cc -Isrc/lib/kissfft -Isrc/lib/kissfft/tools $^ -o $@ -lm

src/test/bin/fptr_cast : src/test/fptr_cast.c
	cc $^ -o $@

src/test/bin/const_types : src/test/const_types.c
	cc $^ -o $@
