source "$(dirname $0)/../scripts/common.sh"
source "$(dirname $0)/attack_examples_common.sh"

[ -z $ATTACK_THRESH ] && ATTACK_THRESH=0.1
[ -z $WINDOW_TYPE ] && WINDOW_TYPE=boxcar
[ -z $ALPHA ] && ALPHA=1
[ -z $H ] && H=256
[ -z $W ] && W=256
export ATTACK_THRESH WINDOW_TYPE ALPHA H W 
INPUT=/tmp/pno2.f64 PYTHONPATH=. python3 test/attack_estimation_spectral_flux.py
