source "$(dirname $0)/../scripts/common.sh"

check_env_set MUMS_PATH "Specify MUMS_PATH"
check_env_set piano_path "Specify piano_path"

[ -z "$SKIP_MUMS_TO_MIDI" ] && \
BASENAME_PREFIX=MPP_PIANO_SOFT_ \
SOURCE_DIR="${MUMS_PATH}/${piano_path}" \
DEST_DIR=/tmp/mums_pno_midi \
bash scripts/mums_to_midi.sh 

[ -z "$SKIP_MUMS_CONV_FMT" ] && \
SOURCE_DIR=/tmp/mums_pno_midi \
DEST_DIR=/tmp/mums_pno_midi_f64 \
bash scripts/mums_conv_fmt.sh 

