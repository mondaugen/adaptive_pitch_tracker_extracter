[ -z $PATH_TO_MUMS ] && PATH_TO_MUMS=~/Documents/sounds/mums/

# yes it's unsafe but it's not a server or anything
tmp_midi_note_name_samples_dir="$(mktemp -u)"

SOURCE_DIR="$PATH_TO_MUMS/disc2/KEYBOARDS/PIANOS/PIANO_MPP_MEDIUM/" \
BASENAME_PREFIX=MPP_PIANO_MEDIUM__ \
DEST_DIR="$tmp_midi_note_name_samples_dir" \
bash scripts/mums_to_midi.sh

SOURCE_DIR="$tmp_midi_note_name_samples_dir" \
DEST_DIR="/tmp/piano" \
bash scripts/mums_conv_fmt.sh

rm -rf "%tmp_midi_note_name_samples_dir"
