.PHONY: do_test build_libdsp_math

R=$(dir $(firstword $(MAKEFILE_LIST)))

do_lookup: $(R)/gen_data.py \
		 $(R)/plot_result.py \
		 $(R)/bin/do_lookup
	PYTHONPATH=. python3 $(R)/gen_data.py && \
	$(R)/bin/do_lookup && \
	PYTHONPATH=. python3 $(R)/plot_result.py

$(R)/bin/do_lookup: $(R)/do_lookup.c \
		src/test/test_common.h \
		build_libdsp_math
	cc -g -Isrc/test -Isrc/lib/dsp_math -Lsrc/lib/dsp_math \
		$(filter %.c,$^) -o $@ -ldsp_math -lm

build_libdsp_math:
	$(MAKE) -C src/lib/dsp_math
