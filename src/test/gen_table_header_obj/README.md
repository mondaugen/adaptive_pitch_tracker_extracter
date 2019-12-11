Generate a table and put directly into object file (without having to
compile), and also make a header so we can access the table and its length.

`src/test/gen_table_header_obj/generate.py` does the generation of binray file and header file
`src/test/gen_table_header_obj/Makefile` wraps resulting binary file to get object file
