#!/bin/bash

# example usage:
# LENGTH=0.2 sox "|scripts/synth_sines.sh %-9 %-7" -d

[ -z $LENGTH ] && LENGTH=2
cmd="sox "
for pitch in $@
do
    cmd+="'"'|sox -n -p synth '"$LENGTH"' sine '"$pitch' "
done
cmd+='-p'
echo "$cmd" | bash
