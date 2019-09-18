# Here's an example of what is possible when time stretching and locking certain
# points

source "$(dirname $0)/../scripts/common.sh"

if [ "$1" == clean ]; then
    rm -rf /tmp/mums_pno_midi /tmp/mums_pno_midi_f64 /tmp/pno.f64 /tmp/pno_stretch_0.1.f64
    exit 0
fi

check_env_set MUMS_PATH "Specify MUMS_PATH"

piano_path=disc2/KEYBOARDS/PIANOS/PIANO_MPP_SOFT/

# prepare the piano samples
bash "$(dirname $0)/tstretch_prepare_piano_samples.sh"

[ -z "$SKIP_SYNTH_AF_SCORE" ] && \
PYTHONPATH=. \
SAMPLE_FILE_DIRECTORY=/tmp/mums_pno_midi_f64 \
SCORE_FILE=test/test_score.txt \
OUTPUT=/tmp/pno.f64 \
python3 scripts/synth_af_score.py 

[ -z "$SKIP_TSTRETCH" ] && \
LOCK_OFFSET=-1024 \
LOCK_WINDOW=2 \
PYTHONPATH=.:scripts \
SCORE_FILE=test/test_score.txt \
SOUND_FILE=/tmp/pno.f64 \
OUTPUT=/tmp/pno_stretch_0.1.f64 \
TSTRETCH_FACTOR=0.1 \
python3 test/tstretch_lock_attacks.py

[ -z "$SKIP_PLAY_ORIGINAL" ] && \
echo "Original" && \
sox -r 16k -c 1 /tmp/pno.f64 -d

[ -z "$SKIP_PLAY_TSTRETCH" ] && \
echo "Time-stretched" && \
sox -r 16k -c 1 /tmp/pno_stretch_0.1.f64 -d
