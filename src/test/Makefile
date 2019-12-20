src/test/bin/kiss_fft_example : src/test/kissfft_example.c \
                                src/lib/kissfft/kiss_fft.c \
                                src/lib/kissfft/tools/kiss_fftr.c
	$(CC) -Isrc/lib/kissfft -Isrc/lib/kissfft/tools $^ -o $@ -lm

src/test/bin/fptr_cast : src/test/fptr_cast.c
	$(CC) $^ -o $@

src/test/bin/const_types : src/test/const_types.c
	$(CC) $^ -o $@

src/test/bin/pvoc_synth_f32_test : src/test/pvoc_synth_f32_test.c \
                                   src/lib/kissfft/kiss_fft.c \
                                   src/lib/kissfft/tools/kiss_fftr.c \
                                   src/lib/datastructures/ola_f32.c \
                                   src/lib/datastructures/windowed_lookup_f32.c \
                                   src/lib/pvoc_synth_f32/routines_linux_native.c \
                                   src/pvoc_synth.c
	$(CC) -g -Isrc/lib/kissfft -Isrc/lib/kissfft/tools -Isrc -Isrc/lib $^ -o $@ -lm

src/test/bin/stretch_f32_file : src/test/stretch_f32_file.c \
                                   src/lib/kissfft/kiss_fft.c \
                                   src/lib/kissfft/tools/kiss_fftr.c \
                                   src/lib/datastructures/ola_f32.c \
                                   src/lib/datastructures/windowed_lookup_f32.c \
                                   src/lib/pvoc_synth_f32/routines_linux_native.c \
                                   src/pvoc_synth.c \
                                   src/test/signal_stretcher_f32.c
	$(CC) -g -Isrc/lib/kissfft -Isrc/lib/kissfft/tools -Isrc -Isrc/lib -Isrc/test $^ -o $@ -lm

src/test/bin/pvoc_synth_ftable_tests_f32 : src/test/pvoc_synth_ftable_tests_f32.c \
                                   src/lib/kissfft/kiss_fft.c \
                                   src/lib/kissfft/tools/kiss_fftr.c \
                                   src/lib/datastructures/ola_f32.c \
                                   src/lib/datastructures/windowed_lookup_f32.c \
                                   src/lib/pvoc_synth_f32/routines_linux_native.c \
                                   src/pvoc_synth.c
	$(CC) -g -Isrc/test -Isrc/lib/kissfft -Isrc/lib/kissfft/tools -Isrc -Isrc/lib $^ -o $@ -lm

src/test/bin/iir_lattice_filter : src/test/iir_lattice_filter.c
	$(CC) -g -Isrc/test $^ -o $@

src/test/bin/iir_lattice_filter_proc : src/test/iir_lattice_filter_proc.c \
                                       src/lib/iir_lattice_filter_f32/routines_linux_native.c
	$(CC) -g -Isrc -Isrc/test $^ -o $@

src/test/bin/lattice_filter_proc : src/test/lattice_filter_proc.c \
                                       src/lib/lattice_filter_f32/routines_linux_native.c
	$(CC) -g -Isrc -Isrc/test $^ -o $@

src/test/bin/adsr_envelopes_test : src/test/adsr_envelopes_test.c \
    src/adsr_envelopes.c \
    src/adsr_envelopes_attack_table/adsr_envelopes_attack_table.o \
    src/lib/adsr_envelopes/adsr_envelopes_linux_native.c \
    src/adsr_envelopes_decay_coeff_table/adsr_envelopes_decay_coeff_table.o
	cc -g -Isrc -Isrc/test -Isrc/adsr_envelopes_attack_table/ \
        -Isrc/adsr_envelopes_decay_coeff_table/ \
        -Isrc/lib/dsp_math/ \
        $^ -o $@ -lm

src/adsr_envelopes_attack_table/adsr_envelopes_attack_table.o : \
    src/adsr_envelopes_attack_table/Makefile \
    src/adsr_envelopes_attack_table/generate_attack_table.py 
	cd src/adsr_envelopes_attack_table/ && $(MAKE)

src/adsr_envelopes_decay_coeff_table/adsr_envelopes_decay_coeff_table.o : \
    src/adsr_envelopes_decay_coeff_table/Makefile \
    src/adsr_envelopes_decay_coeff_table/generate_decay_coeff_table.py
	cd src/adsr_envelopes_decay_coeff_table/ && $(MAKE)

# Basically, profiling such a simple program doesn't work
src/test/bin/ptr_arith_vs_indexing : src/test/profiling/ptr_arith_vs_indexing.o \
    src/test/profiling/aug_note_trans.o
	cc -pg -g $^ -o $@

src/test/profiling/ptr_arith_vs_indexing.o : \
    src/test/profiling/ptr_arith_vs_indexing.c \
    src/test/test_common.h \
    src/test/profiling/aug_note_trans.h
	cc -pg -g -Isrc/test -c $< -o $@

src/test/profiling/aug_note_trans.o : \
    src/test/profiling/aug_note_trans.c \
    src/test/profiling/aug_note_trans.h
	cc -pg -g -O2 -fno-inline -Isrc/test -c $< -o $@

run_note_region_segmenter_test : src/test/bin/note_region_segmenter_test \
    src/test/note_region_segmenter/gen_active_1.py \
    src/test/note_region_segmenter/plot_test.py
	PYTHONPATH=. python3 src/test/note_region_segmenter/gen_active_1.py && \
	$< && \
    PYTHONPATH=. python3 src/test/note_region_segmenter/plot_test.py

run_note_region_segmenter_test_adsr : src/test/bin/note_region_segmenter_test_adsr \
    src/test/note_region_segmenter/gen_active_1.py \
    src/test/note_region_segmenter/plot_test_adsr.py
	PYTHONPATH=. python3 src/test/note_region_segmenter/gen_active_1.py && \
	$< && \
    PYTHONPATH=. python3 src/test/note_region_segmenter/plot_test_adsr.py

src/test/bin/note_region_segmenter_test_adsr : \
    src/test/note_region_segmenter/test_adsr.c \
    src/test/note_region_segmenter.o \
    src/note_region_segmenter.h \
    src/test/test_common.h \
    src/adsr_envelopes.c \
    src/adsr_envelopes_attack_table/adsr_envelopes_attack_table.o \
    src/lib/adsr_envelopes/adsr_envelopes_linux_native.c \
    src/adsr_envelopes_decay_coeff_table/adsr_envelopes_decay_coeff_table.o
	cc -g -Isrc -Isrc/test -Isrc/adsr_envelopes_attack_table/ \
        -Isrc/adsr_envelopes_decay_coeff_table/ \
        -Isrc/lib/dsp_math/ $(filter %.o %.c, $^) -o $@

src/test/bin/note_region_segmenter_test : \
    src/test/note_region_segmenter/test.c \
    src/test/note_region_segmenter.o \
    src/note_region_segmenter.h \
    src/test/test_common.h
	cc -g -Isrc -Isrc/test $(filter %.o %.c, $^) -o $@

src/test/note_region_segmenter.o : src/note_region_segmenter.c \
    src/note_region_segmenter.h
	cc -g -c -Isrc $< -o $@
