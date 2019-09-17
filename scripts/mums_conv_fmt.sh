# Convert a folder of audio files to files of a certain type.
# By default the files are mixed down to mono, converted to a raw f64le file,
# and downsampled to 16KHz sample rate.
# The output mixing and format can be chosen by specifying a different command for CONV.
# The output sample rate can be chosen by specifying SAMPLE_RATE, e.g.,
# SAMPLE_RATE=8000 will resample to 8KHz sample rate.

source "$(dirname $0)/common.sh"

# If you specify a custom CONV statement, you can use \"\$prefix\" to refer to
# the basename of the audio file with the suffix removed, and \"\$DEST_DIR\" the
# output directory, e.g., to specify the filename to output to.
[ -z $CONV ] && CONV=\
'-filter pan='"'"'1c|c0=0.5*c0+0.5*c1'"'"' -f f64le -c:a pcm_f64le \"\$DEST_DIR/\$prefix\".f64'

[ -z $SAMPLE_RATE ] && SAMPLE_RATE=16000

check_env_set SOURCE_DIR "Specify SOURCE_DIR."
check_env_set DEST_DIR "Specify DEST_DIR."

[ ! -d "$DEST_DIR" ] && mkdir -p "$DEST_DIR" && echo "Created directory $DEST_DIR" \
    || fail_msg_exit "Failed creating directory $DEST_DIR"

find "$SOURCE_DIR" -type f | xargs -d'\n' -I{} \
bash -c '
bn=$(basename {})
prefix="${bn%.*}"
eval "ffmpeg -i {} -ar '$SAMPLE_RATE' '"$CONV"'"
'


