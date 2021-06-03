# Tracks using the 'hop' partial tracking method.

# Defaults unless set when calling script
[ -z $PTRACK_DUR ] && export PTRACK_DUR=2.3
[ -z $PTRACK_T0 ] && export PTRACK_T0=0
[ -z $PTRACK_F0 ] && export PTRACK_F0='p:45'
[ -z $SHOW_PLOT ] && export SHOW_PLOT=0
[ -z $PTRACK_OUT ] && export \
    PTRACK_OUT=".testout/sonnet_${PTRACK_T0}_${PTRACK_F0}.f64"

PEAK_T=1e-3 \
PTRACK_REMOVE=0 \
PTRACK_WINTYPE=hann \
PTRACK_METHOD=hop \
PTRACK_WINLEN=2048 \
PTRACK_H=512 \
PTRACK_HARM_LOCK=1 \
PTRACK_F0_MODE=harmonics \
PTRACK_MU=3e-8 \
PTRACK_MAX_STEP=inf \
PTRACK=1 \
PLOT_SD=1 \
SD_BIN_WEIGHT=8 \
SD_H=1024 \
SD_W=4096 \
YSCALE=log \
FMAX=10000 \
TMAX=7.5 \
PYTHONPATH=.:test/gradient_partial_tracker \
NFFT=4096 \
INFILE=sounds/sonnet_m11-12_44.1k_1ch.f64 \
FS=44100 \
SAMPTYPE=float64 \
python3 \
scripts/time_series_specgram.py

