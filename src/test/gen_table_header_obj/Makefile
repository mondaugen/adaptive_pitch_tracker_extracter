.PHONY: clean

use_table.bin : use_table.c table.h table.o
	cc -c $<
	cc $(filter %.o, $^) $(patsubst %.x,%.o,$<) -o $@

table.o : table.f32 create_obj_file.sh
	bash create_obj_file.sh

table.h : generate.py
	python3 $<

table.f32 : generate.py
	python3 $<

clean :
	rm -f table.o table.f32 use_table.bin table.h
