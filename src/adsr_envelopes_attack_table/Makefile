.PHONY: clean

export OUTPUT_FILE_STEM=adsr_envelopes_attack_table
export OUTPUT_TARGET?=elf64-x86-64
export BINARY_ARCHITECTURE?=i386:x86-64

adsr_envelopes_attack_table.o : \
    adsr_envelopes_attack_table.f32 \
    adsr_envelopes_attack_table.h \
    ../scripts/binary_to_obj.sh
	bash ../scripts/binary_to_obj.sh

adsr_envelopes_attack_table.h adsr_envelopes_attack_table.f32 : generate_attack_table.py
	python3 $<

clean :
	rm -f adsr_envelopes_attack_table.h \
        adsr_envelopes_attack_table.f32 \
        adsr_envelopes_attack_table.o

