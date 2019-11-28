src/lib/datastructures/test/bin/ola_f32_test : src/lib/datastructures/test/ola_f32_test.c \
                                               src/lib/datastructures/ola_f32.c
	$(CC) -g -Isrc/lib/datastructures -Isrc/test $^ -o $@

src/lib/datastructures/test/bin/fixed_heap_f32_test : \
    src/lib/datastructures/test/fixed_heap_f32_test.c \
    src/lib/datastructures/fixed_heap.c \
    src/lib/datastructures/fixed_heap_f32.c
	$(CC) -g -I src/test/ -I src/lib/datastructures/ $^ -o $@
