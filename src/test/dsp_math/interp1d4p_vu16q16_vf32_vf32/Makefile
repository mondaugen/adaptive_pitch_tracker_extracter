.PHONY: do_test clean build_libdsp_math

R=$(dir $(firstword $(MAKEFILE_LIST)))

do_test: $(R)/gen_test_data.py \
		 $(R)/chk_ret.py \
		 $(R)/bin/do_interp
	PYTHONPATH=. python3 $(R)/gen_test_data.py && \
	$(R)/bin/do_interp && \
	PYTHONPATH=.:test/interpolation python3 $(R)/chk_ret.py

$(R)/bin/do_interp: $(R)/do_interp.c \
		src/test/test_common.h \
		build_libdsp_math
	cc -g -Isrc/test -Isrc/lib/dsp_math -Lsrc/lib/dsp_math \
		$(filter %.c,$^) -o $@ -ldsp_math -lm

build_libdsp_math:
	$(MAKE) -C src/lib/dsp_math

clean:
	rm -f $(R)/bin/*
