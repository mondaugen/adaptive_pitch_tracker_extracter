dn="$(dirname $0)"
objcopy -I binary -O elf64-x86-64 -B i386:x86-64 --add-symbol table=0 "$dn/table.f32" "$dn/table.o"
