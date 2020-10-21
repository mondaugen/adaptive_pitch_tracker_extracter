.PHONY: runtest

runtest : libadder.so
	LD_LIBRARY_PATH=. python3 adder_so.py

libadder.so : adder.c
	gcc -c -fPIC $< -o $(patsubst %.c, %.o, $<)
	gcc -shared $(patsubst %.c, %.o, $<) -o $@

