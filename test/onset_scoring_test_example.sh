source "$(dirname $0)/../scripts/common.sh"
source "$(dirname $0)/attack_examples_common.sh"

THRESH=.1 \
SOUNDFILE=/tmp/pno2.f64 \
SCOREFILE=test/test_score_2.txt \
PYTHONPATH=.:scripts \
python3 test/onset_scoring_test.py
