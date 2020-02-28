.PHONY: run_closest_index_after_test attacks_from_sd_plot_extraction

src/test/bin/fixed_heap_u32_test : \
    src/test/fixed_heap/u32_test.c \
    src/lib/datastructures/fixed_heap_u32_key.h \
    src/lib/datastructures/fixed_heap_u32_key.c \
    src/lib/datastructures/fixed_heap.h \
    src/lib/datastructures/fixed_heap.c
	cc -g $(filter %.c, $^) -I src/lib/datastructures/ -o $@

run_closest_index_after_test : \
    src/test/bin/closest_index_after_test \
    src/test/attack_finder/closest_index_after_test.py \
    src/test/attack_finder/closest_index_after_test_plot.py
	python3 src/test/attack_finder/closest_index_after_test.py && \
        src/test/bin/closest_index_after_test && \
        PYTHONPATH=. python3 src/test/attack_finder/closest_index_after_test_plot.py

run_closest_index_after_simple_test : \
    src/test/bin/closest_index_after_test \
    src/test/attack_finder/closest_index_after_test.py \
    src/test/attack_finder/closest_index_after_test_plot.py
	python3 src/test/attack_finder/closest_index_after_simple_test.py && \
        src/test/bin/closest_index_after_test && \
        PYTHONPATH=. python3 src/test/attack_finder/closest_index_after_test_plot.py

ATTACK_FINDER_DEPS=src/attack_finder.c \
    src/attack_finder.h \
    src/test/test_common.h \
    src/pvoc_windows.c \
    src/pvoc_windows.h \
    src/find_extrema.c \
    src/find_extrema.h \
    src/lib/datastructures/fixed_heap.h \
    src/lib/datastructures/fixed_heap.c \
    src/lib/datastructures/fixed_heap_u32_key.c \
    src/lib/datastructures/fixed_heap_u32_key.h \
    src/lib/attack_finder/routines_linux_native.c \
    src/lib/dsp_math/dsp_math.h \
    src/lib/dsp_math/dsp_math_linux_native.c \
    src/lib/kissfft/kiss_fft.c \
    src/lib/kissfft/tools/kiss_fftr.c

ATTACK_FINDER_INC=-I src/lib/kissfft/tools/ -I src/lib/kissfft/\
        -Isrc -Isrc/lib/datastructures/ -Isrc/lib/dsp_math 

src/test/bin/closest_index_after_test : \
    src/test/attack_finder/closest_index_after_test.c \
    $(ATTACK_FINDER_DEPS)
	cc -g $(filter %.c, $^) -I src/test $(ATTACK_FINDER_INC) -o $@ -lm

src/test/bin/test_index_mask : \
    src/test/attack_finder/test_index_mask.c \
    $(ATTACK_FINDER_DEPS)
	cc -g $(filter %.c, $^) -I src/test $(ATTACK_FINDER_INC) -o $@ -lm

attacks_from_sd_plot_extraction : \
    src/test/attack_finder/attacks_from_sd_plot_extraction.py \
    src/test/bin/attacks_from_sd
	src/test/bin/attacks_from_sd && \
    PYTHONPATH=. python3 \
    src/test/attack_finder/attacks_from_sd_plot_extraction.py

src/test/bin/attacks_from_sd : \
    src/test/attack_finder/attacks_from_sd_test.c \
    $(ATTACK_FINDER_DEPS)
	cc -g $(filter %.c, $^) -I src/test $(ATTACK_FINDER_INC) -o $@ -lm
