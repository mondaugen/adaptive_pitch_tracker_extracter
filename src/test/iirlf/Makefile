run_adaptive_test :
	make -C ../../lib/iir_lattice_filter_f32/ lib_iirl_f32.so && \
	LD_LIBRARY_PATH=../../lib/iir_lattice_filter_f32/ \
	PYTHONPATH=../../../:../../../test/adaptive_lattice/ \
	python3 ./adaptive_test.py

run_adaptive_test2 :
	make -C ../../lib/iir_lattice_filter_f32/ lib_iirl_f32.so && \
	LD_LIBRARY_PATH=../../lib/iir_lattice_filter_f32/ \
	PYTHONPATH=../../../:../../../test/adaptive_lattice/ \
	python3 ./adaptive_test2.py
