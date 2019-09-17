First you need to prepare some soundfiles. You might want to use MUMS and you might want the MIDI note numbers, in that case you can conver the filenames using

scripts/mums_to_midi.sh

Then they need to be .f64 format, use

scripts/mums_conv_fmt.sh

to convert.

Now you can render a score. An example score can be found at

test/test_score.txt

Here's an example render command line:

PYTHONPATH=~/Documents/development/adaptive_pitch_tracker_extracter OUTPUT=/tmp/render.f64 SCORE_FILE=~/Documents/development/adaptive_pitch_tracker_extracter/test/test_score.txt python3 ~/Documents/development/adaptive_pitch_tracker_extracter/scripts/synth_af_score.py

it assumes you are in the directory so that you needn't specify the path before
the filenames (see the test_score.txt file).
