[ -z $OUTPUT_TARGET ] && OUTPUT_TARGET=elf64-x86-64
[ -z $BINARY_ARCHITECTURE ] && BINARY_ARCHITECTURE=i386:x86-64
[ -z $OUTPUT_FILE_STEM ] && OUTPUT_FILE_STEM=table
dn="$(dirname $0)"
objcopy -I binary -O "$OUTPUT_TARGET" -B "$BINARY_ARCHITECTURE"\
    "$OUTPUT_FILE_STEM.f32" "$OUTPUT_FILE_STEM.o"
objcopy --rename-section .data=.rodata,alloc,load,readonly,data,contents\
    "$OUTPUT_FILE_STEM.o"

