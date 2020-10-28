run_gal_on_file_test :
	make -C ../../src/lib/iir_lattice_filter_f32/ lib_iirl_f32.so && \
	LD_LIBRARY_PATH=../../src/lib/iir_lattice_filter_f32/ \
	PYTHONPATH=../../:. python3 gal_on_file_test.py
