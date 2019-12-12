dn="$(dirname $0)"
#objcopy -I binary -O elf64-x86-64 -B i386:x86-64 --add-symbol table=.rodata:0 "$dn/table.f32" "$dn/table.o"
objcopy -I binary -O elf64-x86-64 -B i386:x86-64 "$dn/table.f32" "$dn/table.o"
objcopy --rename-section .data=.rodata,alloc,load,readonly,data,contents "$dn/table.o"
