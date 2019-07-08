filename=/tmp/sines.f64
sox "|scripts/synth_sines.sh %-9 %-7 %-5 %-4" -t f64 -r 16k -c 1 "$filename"
PYTHONPATH=".:$PYTHONPATH" python3 test/cqt_sines.py

