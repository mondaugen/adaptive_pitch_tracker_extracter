__="$(dirname $0)"

source "$__/common.sh"

python3 "$__/../test/synth_test_midi.py"

synth_midi_to_wav /tmp/midi.mid
