if [ "$1" == clean ]; then
    rm -rf /tmp/mums_pno_midi /tmp/mums_pno_midi_f64 /tmp/pno2.f64
    exit 0
fi

check_env_set MUMS_PATH "Specify MUMS_PATH"

piano_path=disc2/KEYBOARDS/PIANOS/PIANO_MPP_SOFT/

export MUMS_PATH piano_path

# prepare the piano samples
bash "$(dirname $0)/tstretch_prepare_piano_samples.sh"

[ -z "$SKIP_SYNTH_AF_SCORE" ] && \
PYTHONPATH=. \
SAMPLE_FILE_DIRECTORY=/tmp/mums_pno_midi_f64 \
SCORE_FILE=test/test_score_2.txt \
OUTPUT=/tmp/pno2.f64 \
python3 scripts/synth_af_score.py 
