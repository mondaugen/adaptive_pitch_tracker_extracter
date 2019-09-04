# Copy all the files in a folder of MUMS to filenames just containing their MIDI note number
# e.g., ELECTRIC_GUITAR_F#2.mp4 becomes 42.mp4

source "$(dirname $0)/common.sh"

check_env_set SOURCE_DIR "Specify SOURCE_DIR."
check_env_set DEST_DIR "Specify DEST_DIR."
[ -z $BASENAME_PREFIX ] && BASENAME_PREFIX=ELECTRIC_GUITAR_

#find "$SOURCE_DIR" |perl -np -e 's/^ELECTRIC_GUITAR_([A-G]#?\d).*$/$1/'
find "$SOURCE_DIR" -type f | xargs -d'\n' -I{} bash -c \
'
bn=$(basename "{}")
suffix="${bn##*.}"
note=$(echo "$bn"| perl -np -e '"'"'s/^'"$BASENAME_PREFIX"'([A-G]#?\d).*$/$1/'"'"')
midi=$(python3 scripts/note_to_midi.py "$note")
cp "{}" "$DEST_DIR/${midi}.${suffix}"
'
